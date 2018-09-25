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
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================

#ifndef CH_MATERIALSURFACE_SMC_H
#define CH_MATERIALSURFACE_SMC_H

#include "chrono/physics/ChMaterialSurface.h"

namespace chrono {

/// Material data for a surface for use with smooth (penalty) contact method.
/// This data is used to define surface properties owned by ChBody rigid bodies and
/// similar objects; it carries information that is used to make contacts.
class ChApi ChMaterialSurfaceSMC : public ChMaterialSurface {

  public:
    float young_modulus;      ///< Young's modulus (elastic modulus)
    float poisson_ratio;      ///< Poisson ratio
    float static_friction;    ///< Static coefficient of friction
    float sliding_friction;   ///< Kinetic coefficient of friction
    float rolling_friction;   ///< Rolling coefficient of friction
    float twisting_friction;  ///< Twisting coefficient of friction
    float restitution;        ///< Coefficient of restitution
    float constant_adhesion;  ///< Constant adhesion force, when constant adhesion model is used
    float adhesionMultDMT;    ///< Adhesion multiplier used in DMT model.
    float adhesionScheeres;   ///< Adhesion multiplier used in Scheeres model.

    // DMT adhesion model:
    //     adhesion = adhesionMultDMT * sqrt(R_eff).
    // Given the surface energy, w,
    //     adhesionMultDMT = 2 * CH_C_PI * w * sqrt(R_eff).
    // Given the equilibrium penetration distance, y_eq,
    //     adhesionMultDMT = 4.0 / 3.0 * E_eff * powf(y_eq, 1.5)
    // Scheeres adhesion model:
    //     adhesion = adhesionScheeres * r
    // In Scheeres (2010) paper, the following value is used:
    //     adhesion = 3.6 * 10^(-2) S^2 * r
    // with S being the measure of cleanliness

    float kn;  ///< user-specified normal stiffness coefficient
    float kt;  ///< user-specified tangential stiffness coefficient
    float gn;  ///< user-specified normal damping coefficient
    float gt;  ///< user-specified tangential damping coefficient

    ChMaterialSurfaceSMC();
    ChMaterialSurfaceSMC(const ChMaterialSurfaceSMC& other);
    ~ChMaterialSurfaceSMC() {}

    /// "Virtual" copy constructor (covariant return type).
    virtual ChMaterialSurfaceSMC* Clone() const override { return new ChMaterialSurfaceSMC(*this); }

    virtual ContactMethod GetContactMethod() const override { return SMC; }

    /// Young's modulus.
    float GetYoungModulus() const { return young_modulus; }
    void SetYoungModulus(float val) { young_modulus = val; }

    // Poisson ratio.
    float GetPoissonRatio() const { return poisson_ratio; }
    void SetPoissonRatio(float val) { poisson_ratio = val; }

    /// Static and kinetic friction coefficients.
    /// Usually in 0..1 range, rarely above. Default 0.6
    float GetSfriction() const { return static_friction; }
    void SetSfriction(float val) { static_friction = val; }

    float GetKfriction() const { return sliding_friction; }
    void SetKfriction(float val) { sliding_friction = val; }

    /// Rolling friction coefficient.
    /// Usually around 1E-3
    float GetRfriction() const { return rolling_friction; }
    void SetRfriction(float val) { rolling_friction = val; }

	/// Twisting friction coefficient.
	/// Usually around 1E-3
    float GetTfriction() const { return twisting_friction; }
    void SetTfriction(float val) { twisting_friction = val; }
    
    /// Set both static friction and kinetic friction at once, with same value.
    void SetFriction(float val);

    /// Normal restitution coefficient
    float GetRestitution() const { return restitution; }
    void SetRestitution(float val) { 
        double x = val;
        double a0 = 0.81414;
        double a1 = -16.82424;
        double a2 = 134.83788;
        double a3 = -570.11297;
        double a4 = 1520.52153;
        double a5 = -2652.48144;
        double a6 = 3021.53888;
        double a7 = -2165.97278;
        double a8 = 886.75089;
        double a9 = -158.07215;

        double true_restitution;
        true_restitution = a0 + a1 * x + a2 * pow(x,2) + a3 * pow(x,3)
                              + a4 * pow(x,4) + a5 * pow(x,5) + a6 * pow(x,6)
                              + a7 * pow(x,7) + a8 * pow(x,8) + a9 * pow(x,9);
        if (x < 0.1375)
            true_restitution = 0;
        restitution = true_restitution; }

    /// Constant cohesion force
    float GetAdhesion() const { return constant_adhesion; }
    void SetAdhesion(float val) { constant_adhesion = val; }

    /// Adhesion multiplier
    float GetAdhesionMultDMT() const { return adhesionMultDMT; }
    void SetAdhesionMultDMT(float val) { adhesionMultDMT = val; }

    /// Adhesion Scheeres multiplier
    float GetAdhesionScheeres() const { return adhesionScheeres; }
    void SetAdhesionScheeres(float val) { adhesionScheeres = val;}

    /// Stiffness and damping coefficients
    float GetKn() const { return kn; }
    float GetKt() const { return kt; }
    float GetGn() const { return gn; }
    float GetGt() const { return gt; }

    void SetKn(float val) { kn = val; }
    void SetKt(float val) { kt = val; }
    void SetGn(float val) { gn = val; }
    void SetGt(float val) { gt = val; }

    /// Method to allow serializing transient data into in ascii
    /// as a readable item, for example   "chrono::GetLog() << myobject;"
    virtual void StreamOUT(ChStreamOutAscii& mstream) { mstream << "Material SMC \n"; }

    /// Method to allow serialization of transient data to archives.
    virtual void ArchiveOUT(ChArchiveOut& marchive) override;

    /// Method to allow deserialization of transient data from archives.
    virtual void ArchiveIN(ChArchiveIn& marchive) override;
};

CH_CLASS_VERSION(ChMaterialSurfaceSMC, 0)

/// Composite SMC material data for a contact pair.
class ChApi ChMaterialCompositeSMC : public ChMaterialComposite {
  public:
    float E_eff;                ///< Effective elasticity modulus
    float G_eff;                ///< Effective shear modulus
    float mu_eff;               ///< Effective coefficient of friction
    float muR_eff;              ///< Effective coefficient of rolling friction
    float muT_eff;              ///< Effective coefficient of twisting friction
    float cr_eff;               ///< Effective coefficient of restitution
    float adhesion_eff;         ///< Effective cohesion force
    float adhesionMultDMT_eff;  ///< Effective adhesion multiplier (DMT model)
    float adhesionScheeres_eff; ///< Effective adhesion multiplier (Scheeres model)

    float kn;  ///< normal stiffness coefficient
    float kt;  ///< tangential stiffness coefficient
    float gn;  ///< normal viscous damping coefficient
    float gt;  ///< tangential viscuous damping coefficient

    ChMaterialCompositeSMC();

    ChMaterialCompositeSMC(ChMaterialCompositionStrategy<float>* strategy,
                           std::shared_ptr<ChMaterialSurfaceSMC> mat1,
                           std::shared_ptr<ChMaterialSurfaceSMC> mat2);
};

}  // end namespace chrono

#endif
