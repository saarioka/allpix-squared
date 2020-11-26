/**
 * @file
 * @brief perform the parametrisations of carrier mobility depending on material and model
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_MOBILITY_H
#define ALLPIX_MOBILITY_H

#include <stdexcept>
#include <string>
#include <utility>

namespace allpix {

    /**
     * @brief base class for the parametrisation of carrier mobility
     *
     */
    class Mobility {

    public:
        /**
         * @brief constructor
         */
        Mobility() {}

        /**
         * @brief overload the () operator to return the mobility for given carrier type and electric/magnetic
         * field
         */
        double operator()(const CarrierType& type, double efield_mag) { return this->get(type, efield_mag); };

        /**
         * @brief returns value of mobility for a specfic carrier type and electric/magnetic field
         * field
         */
        virtual double get(const CarrierType& type, double efield_mag) = 0;

    protected:
        double electron_Vm_;
        double electron_Ec_;
        double electron_Beta_;

        double hole_Vm_;
        double hole_Ec_;
        double hole_Beta_;

        double electron_A_;
        double electron_Vsat300_;
        double electron_Vsat_;
        double electron_alpha_;
        double electron_p_;
        double electron_mu_zero_;

        double hole_A_;
        double hole_Vsat300_;
        double hole_Vsat_;
        double hole_alpha_;
        double hole_p_;
        double hole_mu_zero_;

        // parameters for compounds
        double x_Ge_;
        double VsA_electron_;
        double VsB_electron_;
        double VsA_hole_;
        double VsB_hole_;
        double VsBow_electron_;
        double VsBow_hole_;
    };

} // namespace allpix

namespace allpix {

    /* @brief class that computes the carrier mobility based on the "Quay" model.
     *
     *"Quay" explained in https://doi.org/10.1016/0038-1101(87)90063-3 uses a parametrization of the
     *Saturation velocity VSat taken from https://doi.org/10.1016/S1369-8001_00)00015-9. This parametrisation has the
     *advantage that it can be applied to most common semi conductor materials.  In this case, the  mobility is a
     *function of VSat and the critical field Ec. Ec is defined as Vsat divided by the mobility at zero field
     *Ec=Vsat/mu_zero. The mobility at zero field is parametrised as mu_zero=alpha*temperature^-p. The parameters alpha
     *(mobility at 0K) and p have to be taken from parametrizations to data. The book Landolt-B�rnstein - Group III
     *Condensed Matterhttps://doi.org/10.1007/b80447 provides a good general summary of the parameters. N.B: the
     *parameter p is generally between 1.5 and 2.3, translating theoretical predictions that mobility at low temperatures
     *(below ~100K) generally goes as T^-1.5 dominated by acoustic phonon scattering, while the mobility at higher
     *temperatures (which goes as T^-2.3) has a contribution from both acoustical and optical phonon scattering.
     */

    class MobilityQuay : public Mobility {

    public:
        MobilityQuay(SensorMaterial detector_material, double temperature) : Mobility() {

            if(detector_material == SensorMaterial::SILICON) {
                electron_A_ = 0.26;
                electron_Vsat300_ = Units::get(1.02e7, "cm/s");
                hole_A_ = 0.37;
                hole_Vsat300_ = Units::get(0.72e7, "cm/s");

                // parameters for mobility at zero field defined in Jacoboni et al
                // https://doi.org/10.1016/0038-1101(77)90054-5
                electron_alpha_ = Units::get(1.43e9, "cm*cm*K/V/s");
                electron_p_ = 2.42;
                hole_alpha_ = Units::get(1.35e8, "cm*cm*K/V/s");
                hole_p_ = 2.20;
                electron_Vsat_ = electron_Vsat300_ / ((1 - electron_A_) + electron_A_ * (temperature / 300));
                hole_Vsat_ = hole_Vsat300_ / ((1 - hole_A_) + hole_A_ * (temperature / 300));
            } else if(detector_material == SensorMaterial::GERMANIUM) {
                electron_A_ = 0.45;
                electron_Vsat300_ = Units::get(0.7e7, "cm/s");
                hole_A_ = 0.39;
                hole_Vsat300_ = Units::get(0.63e7, "cm/s");

                // parameters for mobility at zero field defined in Omar
                // et al for electrons
                // https://doi.org/10.1016/0038-1101(87)90063-3 and in the
                // book Landolt-B�rnstein - Group III Condensed Matter https://doi.org/10.1007/b80447 for holes
                electron_alpha_ = Units::get(5.66e7, "cm*cm*K/V/s");
                electron_p_ = 1.68;
                hole_alpha_ = Units::get(1.05e9, "cm*cm*K/V/s");
                hole_p_ = 2.33;
                electron_Vsat_ = electron_Vsat300_ / ((1 - electron_A_) + electron_A_ * (temperature / 300));

                hole_Vsat_ = hole_Vsat300_ / ((1 - hole_A_) + hole_A_ * (temperature / 300));

            } else if(detector_material == SensorMaterial::SIGE) {
                x_Ge_ = 0.8;                                   // fraction of Silicon hard coded 20%
                VsBow_electron_ = Units::get(-1.28e7, "cm/s"); // The bowing parameter is introduced to account for a
                                                               // nonlinear bowing prevailing in some materials.
                VsBow_hole_ = 0;

                electron_alpha_ = Units::get(4.66e7, "cm*cm*K/V/s");
                electron_p_ = 1.68;

                hole_alpha_ = Units::get(1.05e9, "cm*cm*K/V/s");
                hole_p_ = 2.33;

                VsA_electron_ = Units::get(1.02e7, "cm/s") / ((1 - 0.26) + 0.26 * (temperature / 300.)); // Si electron
                VsB_electron_ = Units::get(0.7e7, "cm/s") / ((1 - 0.45) + 0.45 * (temperature / 300.));  // Ge electron

                VsA_hole_ = Units::get(0.72e7, "cm/s") / ((1 - 0.26) + 0.26 * (temperature / 300.)); // Si hole
                VsB_hole_ = Units::get(0.63e7, "cm/s") / ((1 - 0.45) + 0.45 * (temperature / 300.)); // Ge hole

                electron_Vsat_ =
                    x_Ge_ * VsA_electron_ + (1. - x_Ge_) * VsB_electron_ + x_Ge_ * (1. - x_Ge_) * VsBow_electron_;
                hole_Vsat_ = x_Ge_ * VsA_hole_ + (1 - x_Ge_) * VsB_hole_ + x_Ge_ * (1 - x_Ge_) * VsBow_hole_;
            }

            else if(detector_material == SensorMaterial::GAAS) {
                electron_A_ = 0.44;
                electron_Vsat300_ = Units::get(0.72e7, "cm/s");
                hole_A_ = 0.59;
                hole_Vsat300_ = Units::get(0.9e7, "cm/s");

                // parameters for mobility at zero field defined in Omar
                // et al for electrons
                // https://doi.org/10.1016/0038-1101(87)90063-3 and in the
                // book Landolt-Bornstein - Group III Condensed Matter https://doi.org/10.1007/b80447 for holes
                electron_alpha_ = Units::get(2.5e6, "cm*cm*K/V/s");
                electron_p_ = 1.;
                hole_alpha_ = Units::get(6.3e7, "cm*cm*K/V/s");
                hole_p_ = 2.1;
                electron_Vsat_ = electron_Vsat300_ / ((1 - electron_A_) + electron_A_ * (temperature / 300));
                hole_Vsat_ = hole_Vsat300_ / ((1 - hole_A_) + hole_A_ * (temperature / 300));
            }

            else {
                LOG(INFO) << "Sensor material is \"" << detector_material << "\"";
                throw std::invalid_argument("This mobility parametrization is not valid for the given detector material");
            }
            electron_mu_zero_ = electron_alpha_ / (std::pow(temperature, electron_p_));
            electron_Ec_ = electron_Vsat_ / electron_mu_zero_;

            hole_mu_zero_ = hole_alpha_ / (std::pow(temperature, hole_p_));
            hole_Ec_ = hole_Vsat_ / hole_mu_zero_;
        };

        /**
         * @brief returns the value of mobility for a given electric/magnetic field.
         */
        double get(const CarrierType& type, double efield_mag) {
            double numerator, denominator;
            if(type == CarrierType::ELECTRON) {
                numerator = electron_Vsat_ / electron_Ec_;
                denominator = std::pow(1. + std::pow(efield_mag / electron_Ec_, 2.0), 1.0 / 2.0);
            } else {
                numerator = hole_Vsat_ / hole_Ec_;
                denominator = std::pow(1. + std::pow(efield_mag / hole_Ec_, 2.0), 1.0 / 2.0);
            }
            return numerator / denominator;
        };
    };

    /*@brief class that computes the carrier mobility based on the "Jacoboni" model.
     *The "Jacoboni" model is explained in https://doi.org/10.1016/0038-1101(77)90054-5 (section 5.2). The mobility is
     *a function of a parameter called Vm and the critical field E_c. It is only computed for Silicon.
     */
    class MobilityJacoboni : public Mobility {

    public:
        MobilityJacoboni(SensorMaterial detector_material, double temperature) : Mobility() {
            if(detector_material == SensorMaterial::SILICON) // this parametrisation is only available for Silicon
            {
                electron_Vm_ = Units::get(1.53e9 * std::pow(temperature, -0.87), "cm/s");
                electron_Ec_ = Units::get(1.01 * std::pow(temperature, 1.55), "V/cm");
                electron_Beta_ = 2.57e-2 * std::pow(temperature, 0.66);
                hole_Vm_ = Units::get(1.62e8 * std::pow(temperature, -0.52), "cm/s");
                hole_Ec_ = Units::get(1.24 * std::pow(temperature, 1.68), "V/cm");
                hole_Beta_ = 0.46 * std::pow(temperature, 0.17);

            } else {
                LOG(INFO) << "Sensor material is \"" << detector_material << "\"";
                throw std::invalid_argument("This mobility parametrization is not valid for the given detector material");
            }
        };

        /**
         * @brief returns the value of mobility for a given electric/magnetic field.
         */
        double get(const CarrierType& type, double efield_mag) {
            double numerator, denominator;
            // Parameterization variables for Silicon from https://doi.org/10.1016/0038-1101(77)90054-5 (section 5.2)
            if(type == CarrierType::ELECTRON) {
                numerator = electron_Vm_ / electron_Ec_;
                denominator = std::pow(1. + std::pow(efield_mag / electron_Ec_, electron_Beta_), 1.0 / electron_Beta_);
            } else {
                numerator = hole_Vm_ / hole_Ec_;
                denominator = std::pow(1. + std::pow(efield_mag / hole_Ec_, hole_Beta_), 1.0 / hole_Beta_);
            }
            return numerator / denominator;
        };
    };
} // namespace allpix

#endif
