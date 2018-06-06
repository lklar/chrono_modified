// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban, Alessandro Tasora
// =============================================================================
//
// Smooth (penalty-based) contact between two generic contactable objects.
//
// =============================================================================

#ifndef CH_CONTACT_SMC_H
#define CH_CONTACT_SMC_H

#include <algorithm>
#include <cmath>

#include "chrono/collision/ChCCollisionModel.h"
#include "chrono/core/ChFrame.h"
#include "chrono/core/ChMatrixDynamic.h"
#include "chrono/core/ChVectorDynamic.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChContactTuple.h"
#include "chrono/physics/ChMaterialSurfaceSMC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChKblockGeneric.h"
#include "chrono/solver/ChSystemDescriptor.h"
#include "chrono/timestepper/ChState.h"

namespace chrono {

/// Class for smooth (penalty-based) contact between two generic contactable objects.
/// Ta and Tb are of ChContactable sub classes.
template <class Ta, class Tb>
class ChContactSMC : public ChContactTuple<Ta, Tb> {
  public:
    typedef typename ChContactTuple<Ta, Tb>::typecarr_a typecarr_a;
    typedef typename ChContactTuple<Ta, Tb>::typecarr_b typecarr_b;

  private:
    struct ChContactJacobian {
        ChKblockGeneric m_KRM;        ///< sum of scaled K and R, with pointers to sparse variables
        ChMatrixDynamic<double> m_K;  ///< K = dQ/dx
        ChMatrixDynamic<double> m_R;  ///< R = dQ/dv
    };

    ChVector<> m_force;        ///< contact force on objB
    ChContactJacobian* m_Jac;  ///< contact Jacobian data
    double kn, gn, kt, gt, effRadius;

  public:
    ChContactSMC() : m_Jac(NULL) {}

    ChContactSMC(ChContactContainer* mcontainer,          ///< contact container
                 Ta* mobjA,                               ///< collidable object A
                 Tb* mobjB,                               ///< collidable object B
                 const collision::ChCollisionInfo& cinfo  ///< data for the contact pair
                 )
        : ChContactTuple<Ta, Tb>(mcontainer, mobjA, mobjB, cinfo), m_Jac(NULL) {
        SetParameters(mobjA, mobjB, cinfo);
        Reset(mobjA, mobjB, cinfo);
    }

    ~ChContactSMC() { delete m_Jac; }

    /// Get the contact force, if computed, in contact coordinate system
    virtual ChVector<> GetContactForce() const override { return this->contact_plane.MatrT_x_Vect(m_force); }

    /// Get the contact penetration (positive if there is overlap).
    double GetContactPenetration() const { return -this->norm_dist; }

    /// Get the contact force, expressed in the frame of the contact.
    ChVector<> GetContactForceAbs() const { return m_force; }

    /// Access the proxy to the Jacobian.
    const ChKblockGeneric* GetJacobianKRM() const { return m_Jac ? &(m_Jac->m_KRM) : NULL; }
    const ChMatrixDynamic<double>* GetJacobianK() const { return m_Jac ? &(m_Jac->m_K) : NULL; }
    const ChMatrixDynamic<double>* GetJacobianR() const { return m_Jac ? &(m_Jac->m_R) : NULL; }

    /// Reinitialize this contact.
    virtual void Reset(Ta* mobjA,                               ///< collidable object A
                       Tb* mobjB,                               ///< collidable object B
                       const collision::ChCollisionInfo& cinfo  ///< data for the contact pair
                       ) override {
        // Inherit base class.
        ChContactTuple<Ta, Tb>::Reset(mobjA, mobjB, cinfo);

        // Note: cinfo.distance is the same as this->norm_dist.
        assert(cinfo.distance < 0);

        // Calculate composite material properties
        ChMaterialCompositeSMC mat(
            this->container->GetSystem()->composition_strategy.get(),
            std::static_pointer_cast<ChMaterialSurfaceSMC>(this->objA->GetMaterialSurfaceBase()),
            std::static_pointer_cast<ChMaterialSurfaceSMC>(this->objB->GetMaterialSurfaceBase()));

        // Check for a user-provided callback to modify the material
        if (this->container->GetAddContactCallback()) {
            this->container->GetAddContactCallback()->OnAddContact(cinfo, &mat);
        }

        // Calculate contact force.
        m_force = CalculateForce(-this->norm_dist,                            // overlap (here, always positive)
                                 this->normal,                                // normal contact direction
                                 this->objA->GetContactPointSpeed(this->p1),  // velocity of contact point on objA
                                 this->objB->GetContactPointSpeed(this->p2),  // velocity of contact point on objB
                                 mat                                          // composite material for contact pair
        );

        // Set up and compute Jacobian matrices.
        if (static_cast<ChSystemSMC*>(this->container->GetSystem())->GetStiffContact()) {
            CreateJacobians();
            CalculateJacobians(mat);
        }
    }

    /// Calculate coefficents of stiffness and dampening, effective radius of curvature and pre-contact velocity for
    /// this contact
    void SetParameters(Ta* mobjA, Tb* mobjB, const collision::ChCollisionInfo& cinfo) {
        // Calculate effective radius of curvature by deducing the radius of the objects (assuming theyre spheres)
        // from the dimensions of their collision model
        ChVector<> bbminA, bbmaxA, bbminB, bbmaxB;
        cinfo.modelA->GetAABB(bbminA, bbmaxA);
        cinfo.modelB->GetAABB(bbminB, bbmaxB);
        double radA = ((bbmaxA - bbminA).x() / 2.0) - cinfo.modelA->GetEnvelope(),
               radB = ((bbmaxB - bbminB).x() / 2.0) - cinfo.modelB->GetEnvelope();
        effRadius = radA * radB / (radA + radB);

        // Calculate normal relative velocity at the start of the contact
        ChVector<> velocity1 = this->objA->GetContactPointSpeed(this->p1),
                   velocity2 = this->objB->GetContactPointSpeed(this->p2), normalDir = this->normal;
        ChVector<> relvel = velocity2 - velocity1;
        double pre_v_rel_mag = abs(relvel.Dot(normalDir));

        // Extract parameters from containing system
        ChSystemSMC* sys = static_cast<ChSystemSMC*>(this->container->GetSystem());
        bool use_mat_props = sys->UsingMaterialProperties();
        ChSystemSMC::ContactForceModel contact_model = sys->GetContactForceModel();

        // Calculate effective mass
        double eff_mass = this->objA->GetContactableMass() * this->objB->GetContactableMass() /
                          (this->objA->GetContactableMass() + this->objB->GetContactableMass());

        // Calculate composite material properties
        ChMaterialCompositeSMC mat(
            this->container->GetSystem()->composition_strategy.get(),
            std::static_pointer_cast<ChMaterialSurfaceSMC>(this->objA->GetMaterialSurfaceBase()),
            std::static_pointer_cast<ChMaterialSurfaceSMC>(this->objB->GetMaterialSurfaceBase()));

        // Check for a user-provided callback to modify the material
        if (this->container->GetAddContactCallback()) {
            this->container->GetAddContactCallback()->OnAddContact(cinfo, &mat);
        }

        double delta = -this->norm_dist;

        // Calculate stiffness and viscous damping coefficients.

        switch (contact_model) {
            case ChSystemSMC::Hooke:
                if (use_mat_props) {
                    double tmp_k = (16.0 / 15) * std::sqrt(effRadius) * mat.E_eff;
                    double v2 = sys->GetCharacteristicImpactVelocity() * sys->GetCharacteristicImpactVelocity();
                    double loge = (mat.cr_eff < CH_MICROTOL) ? std::log(CH_MICROTOL) : std::log(mat.cr_eff);
                    loge = (mat.cr_eff > 1 - CH_MICROTOL) ? std::log(1 - CH_MICROTOL) : loge;
                    double tmp_g = 1 + std::pow(CH_C_PI / loge, 2);
                    kn = tmp_k * std::pow(eff_mass * v2 / tmp_k, 1.0 / 5);
                    kt = kn;
                    gn = std::sqrt(4 * eff_mass * kn / tmp_g);
                    gt = gn;
                } else {
                    kn = mat.kn;
                    kt = mat.kt;
                    gn = eff_mass * mat.gn;
                    gt = eff_mass * mat.gt;
                }

                break;

            case ChSystemSMC::Hertz:
                if (use_mat_props) {
                    double sqrt_R = std::sqrt(effRadius);
                    double Sn = 2 * mat.E_eff * sqrt_R;
                    double St = 8 * mat.G_eff * sqrt_R;
                    double loge = (mat.cr_eff < CH_MICROTOL) ? std::log(CH_MICROTOL) : std::log(mat.cr_eff);
                    double beta = loge / std::sqrt(loge * loge + CH_C_PI * CH_C_PI);
                    kn = (2.0 / 3) * Sn;
                    kt = St;
                    gn = -2 * std::sqrt(5.0 / 6) * beta * std::sqrt(Sn * eff_mass);
                    gt = -2 * std::sqrt(5.0 / 6) * beta * std::sqrt(St * eff_mass);
                } else {
                    double tmp = effRadius * std::sqrt(delta);
                    kn = tmp * mat.kn;
                    kt = tmp * mat.kt;
                    gn = tmp * eff_mass * mat.gn;
                    gt = tmp * eff_mass * mat.gt;
                }

                break;

            case ChSystemSMC::Flores:
                if (use_mat_props) {
                    double cor = (mat.cr_eff < CH_MICROTOL) ? CH_MICROTOL : mat.cr_eff;
                    cor = (mat.cr_eff > 1 - CH_MICROTOL) ? 1 - CH_MICROTOL : cor;
                    kn = (4.0 / 3.0) * mat.E_eff * std::sqrt(effRadius);
                    kt = kn;
                    gn = 8.0 * (1.0 - cor) / (5.0 * cor * pre_v_rel_mag);
                    gt = gn;
                } else {
                    double tmp = eff_radius * std::sqrt(delta);
                    kn = tmp * mat.kn;
                    kt = tmp * mat.kt;
                    gn = tmp * eff_mass * mat.gn;
                    gt = tmp * eff_mass * mat.gt;
                }

                break;

            case ChSystemSMC::PlainCoulomb:
                if (use_mat_props) {
                    double Sn = 2 * mat.E_eff;
                    double St = 8 * mat.G_eff;
                    double loge = (mat.cr_eff < CH_MICROTOL) ? std::log(CH_MICROTOL) : std::log(mat.cr_eff);
                    double beta = loge / std::sqrt(loge * loge + CH_C_PI * CH_C_PI);
                    kn = (2.0 / 3) * Sn;
                    gn = -2 * std::sqrt(5.0 / 6) * beta * std::sqrt(Sn * eff_mass);
                } else {
                    double tmp = std::sqrt(delta);
                    kn = tmp * mat.kn;
                    gn = tmp * mat.gn;
                }

                kt = 0;
                gt = 0;
        }
    }

    /// Calculate contact force, expressed in absolute coordinates.
    ChVector<> CalculateForce(
        double delta,                      ///< overlap in normal direction
        const ChVector<>& normal_dir,      ///< normal contact direction (expressed in global frame)
        const ChVector<>& vel1,            ///< velocity of contact point on objA (expressed in global frame)
        const ChVector<>& vel2,            ///< velocity of contact point on objB (expressed in global frame)
        const ChMaterialCompositeSMC& mat  ///< composite material for contact pair
    ) {
        // Set contact force to zero if no penetration.
        if (delta <= 0) {
            return ChVector<>(0, 0, 0);
        }

        // Extract parameters from containing system
        ChSystemSMC* sys = static_cast<ChSystemSMC*>(this->container->GetSystem());
        double dT = sys->GetStep();
        bool use_mat_props = sys->UsingMaterialProperties();
        ChSystemSMC::ContactForceModel contact_model = sys->GetContactForceModel();
        ChSystemSMC::AdhesionForceModel adhesion_model = sys->GetAdhesionForceModel();
        ChSystemSMC::TangentialDisplacementModel tdispl_model = sys->GetTangentialDisplacementModel();

        // Relative velocity at contact
        ChVector<> relvel = vel2 - vel1;
        double relvel_n_mag = relvel.Dot(normal_dir);
        ChVector<> relvel_n = relvel_n_mag * normal_dir;
        ChVector<> relvel_t = relvel - relvel_n;
        double relvel_t_mag = relvel_t.Length();

        // Tangential displacement (magnitude)
        double delta_t = 0;
        switch (tdispl_model) {
            case ChSystemSMC::OneStep:
                delta_t = relvel_t_mag * dT;
                break;
            case ChSystemSMC::MultiStep:
                //// TODO: implement proper MultiStep mode
                delta_t = relvel_t_mag * dT;
                break;
            default:
                break;
        }

        double forceN, forceT;

        // Calculate the magnitudes of the normal and tangential contact forces
        switch (contact_model) {
            case ChSystemSMC::Hooke:
           // case ChSystemSMC::Flores:
                forceN = kn * delta - gn * relvel_n_mag;
                forceT = kt * delta_t - gt * relvel_t_mag;
                break;

			case ChSystemSMC::Hertz:
            case ChSystemSMC::PlainCoulomb:
                forceN = kn * pow(delta, 3.0 / 2.0) - gn * pow(delta, 1.0 / 4.0) * relvel_n_mag;
                forceT = kt * pow(delta_t, 3.0 / 2.0) - gt * pow(delta_t, 1.0 / 4.0) * relvel_t_mag;
                break;

            case ChSystemSMC::Flores:
                forceN = kn * pow(delta, 3.0 / 2.0) * (1 - gn * relvel_n_mag);
                forceT = kt * pow(delta, 3.0 / 2.0) * (1 - gn * relvel_n_mag);
                break;
        }

        // If the resulting normal contact force is negative, the two shapes are moving
        // away from each other so fast that no contact force is generated.
        if (forceN < 0) {
            forceN = 0;
            forceT = 0;
        }

        // Include adhesion force
        switch (adhesion_model) {
            case ChSystemSMC::Constant:
                forceN -= mat.adhesion_eff;
                break;
            case ChSystemSMC::DMT:
                forceN -= mat.adhesionMultDMT_eff * sqrt(effRadius);
                break;
        }

        // Coulomb law
        forceT = std::min<double>(forceT, mat.mu_eff * std::abs(forceN));

        // Accumulate normal and tangential forces
        ChVector<> force = forceN * normal_dir;
        if (relvel_t_mag >= sys->GetSlipVelocitythreshold())
            force -= (forceT / relvel_t_mag) * relvel_t;

        return force;
    }

    /// Compute all forces in a contiguous array.
    /// Used in finite-difference Jacobian approximation.
    void CalculateQ(const ChState& stateA_x,            ///< state positions for objA
                    const ChStateDelta& stateA_w,       ///< state velocities for objA
                    const ChState& stateB_x,            ///< state positions for objB
                    const ChStateDelta& stateB_w,       ///< state velocities for objB
                    const ChMaterialCompositeSMC& mat,  ///< composite material for contact pair
                    ChVectorDynamic<>& Q                ///< output generalized forces
    ) {
        // Express contact points in local frames.
        // We assume that these points remain fixed to their respective contactable objects.
        ChVector<> p1_loc = this->objA->GetCsysForCollisionModel().TransformPointParentToLocal(this->p1);
        ChVector<> p2_loc = this->objB->GetCsysForCollisionModel().TransformPointParentToLocal(this->p2);

        // Express the local points in global frame
        ChVector<> p1_abs = this->objA->GetContactPoint(p1_loc, stateA_x);
        ChVector<> p2_abs = this->objB->GetContactPoint(p2_loc, stateB_x);

        /*
            Note: while this can be somewhat justified for a ChBody, it will not work
                  for a mesh vertex for instance...

        // Project the points onto the unperturbed normal line
        p1_abs = this->p1 + Vdot(p1_abs - this->p1, this->normal) * this->normal;
        p2_abs = this->p2 + Vdot(p2_abs - this->p2, this->normal) * this->normal;
        */

        // Calculate normal direction (expressed in global frame)
        ChVector<> normal_dir = (p1_abs - p2_abs).GetNormalized();

        // Calculate penetration depth
        double delta = (p1_abs - p2_abs).Length();

        // If the normal direction flipped sign, change sign of delta
        if (Vdot(normal_dir, this->normal) < 0)
            delta = -delta;

        // Calculate velocity of contact points (expressed in global frame)
        ChVector<> vel1 = this->objA->GetContactPointSpeed(p1_loc, stateA_x, stateA_w);
        ChVector<> vel2 = this->objB->GetContactPointSpeed(p2_loc, stateB_x, stateB_w);

        // Compute the contact force.
        ChVector<> force = CalculateForce(delta, normal_dir, vel1, vel2, mat);

        // Compute and load the generalized contact forces.
        this->objA->ContactForceLoadQ(-force, p1_abs, stateA_x, Q, 0);
        this->objB->ContactForceLoadQ(force, p2_abs, stateB_x, Q, this->objA->ContactableGet_ndof_w());
    }

    /// Create the Jacobian matrices.
    /// These matrices are created/resized as needed.
    void CreateJacobians() {
        delete m_Jac;
        m_Jac = new ChContactJacobian;

        // Set variables and resize Jacobian matrices.
        // NOTE: currently, only contactable objects derived from ChContactable_1vars<6>,
        //       ChContactable_1vars<3>, and ChContactable_3vars<3,3,3> are supported.
        int ndof_w = 0;
        std::vector<ChVariables*> vars;

        vars.push_back(this->objA->GetVariables1());
        if (auto objA_333 = dynamic_cast<ChContactable_3vars<3, 3, 3>*>(this->objA)) {
            vars.push_back(objA_333->GetVariables2());
            vars.push_back(objA_333->GetVariables3());
        }
        ndof_w += this->objA->ContactableGet_ndof_w();

        vars.push_back(this->objB->GetVariables1());
        if (auto objB_333 = dynamic_cast<ChContactable_3vars<3, 3, 3>*>(this->objB)) {
            vars.push_back(objB_333->GetVariables2());
            vars.push_back(objB_333->GetVariables3());
        }
        ndof_w += this->objB->ContactableGet_ndof_w();

        m_Jac->m_KRM.SetVariables(vars);
        m_Jac->m_K.Reset(ndof_w, ndof_w);
        m_Jac->m_R.Reset(ndof_w, ndof_w);
        assert(m_Jac->m_KRM.Get_K()->GetColumns() == ndof_w);
    }

    /// Calculate Jacobian of generalized contact forces.
    void CalculateJacobians(const ChMaterialCompositeSMC& mat) {
        // Compute a finite-difference approximations to the Jacobians of the contact forces and
        // load dQ/dx into m_Jac->m_K and dQ/dw into m_Jac->m_R.
        // Note that we only calculate these Jacobians whenever the contact force itself is calculated,
        // that is only once per step.  The Jacobian of generalized contact forces will therefore be
        // constant over the time step.

        // Get states for objA
        int ndofA_x = this->objA->ContactableGet_ndof_x();
        int ndofA_w = this->objA->ContactableGet_ndof_w();
        ChState stateA_x(ndofA_x, NULL);
        ChStateDelta stateA_w(ndofA_w, NULL);
        this->objA->ContactableGetStateBlock_x(stateA_x);
        this->objA->ContactableGetStateBlock_w(stateA_w);

        // Get states for objB
        int ndofB_x = this->objB->ContactableGet_ndof_x();
        int ndofB_w = this->objB->ContactableGet_ndof_w();
        ChState stateB_x(ndofB_x, NULL);
        ChStateDelta stateB_w(ndofB_w, NULL);
        this->objB->ContactableGetStateBlock_x(stateB_x);
        this->objB->ContactableGetStateBlock_w(stateB_w);

        // Compute Q at current state
        ChVectorDynamic<> Q0(ndofA_w + ndofB_w);
        CalculateQ(stateA_x, stateA_w, stateB_x, stateB_w, mat, Q0);

        // Finite-difference approximation perturbation.
        // Note that ChState and ChStateDelta are set to 0 on construction.
        // To accommodate objects with quaternion states, use the method ContactableIncrementState while
        // calculating Jacobian columns corresponding to position states.
        double perturbation = 1e-5;
        ChState stateA_x1(ndofA_x, NULL);
        ChState stateB_x1(ndofB_x, NULL);
        ChStateDelta prtrbA(ndofA_w, NULL);
        ChStateDelta prtrbB(ndofB_w, NULL);

        ChVectorDynamic<> Q1(ndofA_w + ndofB_w);
        ChVectorDynamic<> Jcolumn(ndofA_w + ndofB_w);

        // Jacobian w.r.t. variables of objA
        for (int i = 0; i < ndofA_w; i++) {
            prtrbA(i) += perturbation;
            this->objA->ContactableIncrementState(stateA_x, prtrbA, stateA_x1);
            CalculateQ(stateA_x1, stateA_w, stateB_x, stateB_w, mat, Q1);
            prtrbA(i) -= perturbation;

            Jcolumn = (Q1 - Q0) * (-1 / perturbation);  // note sign change
            m_Jac->m_K.PasteMatrix(Jcolumn, 0, i);

            stateA_w(i) += perturbation;
            CalculateQ(stateA_x, stateA_w, stateB_x, stateB_w, mat, Q1);
            stateA_w(i) -= perturbation;

            Jcolumn = (Q1 - Q0) * (-1 / perturbation);  // note sign change
            m_Jac->m_R.PasteMatrix(Jcolumn, 0, i);
        }

        // Jacobian w.r.t. variables of objB
        for (int i = 0; i < ndofB_w; i++) {
            prtrbB(i) += perturbation;
            this->objB->ContactableIncrementState(stateB_x, prtrbB, stateB_x1);
            CalculateQ(stateA_x, stateA_w, stateB_x1, stateB_w, mat, Q1);
            prtrbB(i) -= perturbation;

            Jcolumn = (Q1 - Q0) * (-1 / perturbation);  // note sign change
            m_Jac->m_K.PasteMatrix(Jcolumn, 0, ndofA_w + i);

            stateB_w(i) += perturbation;
            CalculateQ(stateA_x, stateA_w, stateB_x, stateB_w, mat, Q1);
            stateB_w(i) -= perturbation;

            Jcolumn = (Q1 - Q0) * (-1 / perturbation);  // note sign change
            m_Jac->m_R.PasteMatrix(Jcolumn, 0, ndofA_w + i);
        }
    }

    /// Apply contact forces to the two objects.
    /// (new version, for interfacing to ChTimestepper and ChIntegrable)
    virtual void ContIntLoadResidual_F(ChVectorDynamic<>& R, const double c) override {
        ChVector<> abs_force_scaled(m_force * c);

        if (this->objA->IsContactActive())
            this->objA->ContactForceLoadResidual_F(-abs_force_scaled, this->p1, R);

        if (this->objB->IsContactActive())
            this->objB->ContactForceLoadResidual_F(abs_force_scaled, this->p2, R);
    }

    /// Inject Jacobian blocks into the system descriptor.
    /// Tell to a system descriptor that there are item(s) of type ChKblock in this object
    /// (for further passing it to a solver)
    virtual void ContInjectKRMmatrices(ChSystemDescriptor& mdescriptor) override {
        if (m_Jac)
            mdescriptor.InsertKblock(&m_Jac->m_KRM);
    }

    /// Compute Jacobian of contact forces.
    virtual void ContKRMmatricesLoad(double Kfactor, double Rfactor) override {
        if (m_Jac) {
            m_Jac->m_KRM.Get_K()->FillElem(0);

            m_Jac->m_KRM.Get_K()->MatrInc(m_Jac->m_K * Kfactor);
            m_Jac->m_KRM.Get_K()->MatrInc(m_Jac->m_R * Rfactor);
        }
    }
};

}  // end namespace chrono

#endif
