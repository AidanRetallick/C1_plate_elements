//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Header file for the Biharmonic Bell elements
#ifndef OOMPH_MYBIHARMONIC_ELEMENTS_HEADER
#define OOMPH_MYBIHARMONIC_ELEMENTS_HEADER


#include<sstream>

//OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"


namespace oomph
{
//=============================================================
/// A class for all subparametric elements that solve the 2D-
/// Biharmonic equations.
/// \f[
/// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
/// \f]
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
class FoepplVonKarmanEquations : public virtual FiniteElement
{

public:

 /// \short Function pointer to an alphaDT type swelling function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*SwellingFctPt)(const Vector<double>& x, double& f);

 /// \short Function pointer to pressure function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*PressureFctPt)(const Vector<double>& x, double& f);


 /// \short Function pointer to in plane forcing function  fct(x,g(x)) --
 /// x is a Vector!
 typedef void (*InPlaneForcingFctPt)(const Vector<double>& x,
                                            Vector<double>& forcing);

 /// \short Function pointer to the Error Metric we are using
 ///  e.g could be that we are just interested in error on w etc.
 typedef void (*ErrorMetricFctPt)(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm);

 /// \short Function pointer to the Error Metric we are using if we want multiple
 ///  errors.  e.g could be we want errors seperately on each displacment
 typedef void (*MultipleErrorMetricFctPt)(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, Vector<double>& error, 
  Vector<double>& norm);

 /// \short (pure virtual) interface to fixing the out of plane displacement
 virtual void fix_out_of_plane_displacement_dof(const unsigned& dof_number, 
const unsigned& boundary_number, const PressureFctPt& w)=0;
 // NB do we need this?

 /// \short (pure virtual) interface to fixing the in plane displacement
 virtual void fix_in_plane_displacement_dof(const unsigned& dof_number,
const unsigned& boundary_number, const PressureFctPt& u)=0;
 // NB do we need this?
 
 /// \short (Read only) Access to number of internal dofs
// unsigned number_of_internal_dofs() const {return this->nbubble_basis();}
//
//
 // Get the number of nodes of the out-of-plane functions, pure virtual
 virtual unsigned nnode_outofplane() const =  0;

 // Get the number of nodes of the out-of-plane functions, pure virtual
 virtual unsigned nnode_inplane() const =  0;

 // Get the number of basis functions, pure virtual
 virtual unsigned nnodal_basis_type() const =  0;

 // Get the number of internal basis functions, pure virtual
 virtual unsigned nbubble_basis() const = 0;

 // Get the number of internal basis functions, pure virtual
 virtual unsigned nbubble_basis_type() const = 0;
 
 public:

 /// Get pointer to association matrix 
 DenseMatrix<double> *get_association_matrix_pt()const {return Association_matrix_pt;};

 public: 

 /// Eta
 const double &eta() const {return *Eta_pt;}

 /// Pointer to eta
 const double* &eta_pt() {return Eta_pt;}

 /// Pure virtual function to pin all deflection dofs
 virtual void pin_all_deflection_dofs() const=0;

 /// Constructor (must initialise the Pressure_fct_pt to null)
 FoepplVonKarmanEquations() : Pressure_fct_pt(0),
   In_plane_forcing_fct_pt(0), Error_metric_fct_pt(0),
    Multiple_error_metric_fct_pt(0), Association_matrix_pt(0) 
  {
   Eta_pt = &Default_Eta_Value;
   Nu_pt = &Default_Nu_Value;
  }

 /// Broken copy constructor
 FoepplVonKarmanEquations(const FoepplVonKarmanEquations& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanEquations");
  }

 /// Broken assignment operator
 void operator=(const FoepplVonKarmanEquations&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanEquations");
  }

 /// \short Return the index at which the first u unknown value
 /// is stored.
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 /// Note these are stored before w_dofs otherwise at vertex nodes the 
 /// stored value would be at a different index
 virtual inline unsigned u_nodal_index_foeppl_von_karman() const {return 0;}

 /// \short Return the index at which the first w unknown value
 /// is stored.
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 /// Note that at midside nodes there are no w dofs 
 virtual inline unsigned w_nodal_index_foeppl_von_karman() const {return u_nodal_index_foeppl_von_karman()+2;}

 /// Output with default number of plot points
 void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   FoepplVonKarmanEquations::output(outfile,n_plot);
  }

 /// \short Output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   FoepplVonKarmanEquations::output(file_pt,n_plot);
  }

 /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);

 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

 /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
 /// n_plot^DIM plot points (dummy time-dependent version to
 /// keep intel compiler happy)
 virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
                         const double& time,
                         FiniteElement::UnsteadyExactSolutionFctPt
                         exact_soln_pt)
  {
   throw OomphLibError(
    "There is no time-dependent output_fct() for these elements ",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm);

 /// Get error against and norm of exact solution
 void compute_error_in_deflection(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm);

 /// Dummy, time dependent error checker
 void compute_error(std::ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)
 {
  throw OomphLibError(
		      "There is no time-dependent compute_error() for these elements",
		      OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
 }

 /// Fill in the strain tensor from displacement gradients
 void get_epsilon(DenseMatrix<double>& epsilon,
		  const DenseMatrix<double>& grad_u,
		  const DenseMatrix<double>& grad_w,
		  const double& c_swell)const
 {
  // Truncated Green Lagrange strain tensor
  DenseMatrix<double> dummy_epsilon(this->dim(),this->dim(),0.0);
  for(unsigned alpha=0;alpha<this->dim();++alpha)
   {
    for(unsigned beta=0;beta<this->dim();++beta)
     {
      // Truncated Green Lagrange strain tensor
      dummy_epsilon(alpha,beta) += 0.5* grad_u(alpha,beta)
       + 0.5*grad_u(beta,alpha)
       + 0.5*grad_w(0,alpha)*grad_w(0,beta);
     }
    // Swelling slack
    dummy_epsilon(alpha,alpha) -= c_swell;
    
   }
  epsilon=dummy_epsilon;
 }

 /// Fill in the stress tensor from displacement gradients
 void get_sigma(DenseMatrix<double>& sigma,
		const DenseMatrix<double>& grad_u,
		const DenseMatrix<double>& grad_w,
		const double c_swell)const
 {
  // Get the Poisson ratio
  double nu(get_nu());

  // Truncated Green Lagrange strain tensor
  DenseMatrix<double> epsilon(this->dim(),this->dim(),0.0);
  get_epsilon(epsilon, grad_u, grad_w, c_swell);
   
  // Now construct the Stress
  for(unsigned alpha=0;alpha<this->dim();++alpha)
   {
    for(unsigned beta=0;beta<this->dim();++beta)
     {
      // The Laplacian term: Trace[ \epsilon ] I
      // \nu * \epsilon_{\alpha \beta} delta_{\gamma \gamma}
      sigma(alpha,alpha) += nu*epsilon(beta,beta)/(1-nu*nu);
      
      // The scalar transform term: \epsilon
      // (1-\nu) * \epsilon_{\alpha \beta}
      sigma(alpha,beta) += (1-nu)* epsilon(alpha,beta)/(1-nu*nu);
     }
   }
 }

 // TODO: Fix this function to use index loops
 /// Fill in the stress tensor using a precalculated strain tensor
 void get_sigma_from_epsilon(DenseMatrix<double>& sigma,
			     const DenseMatrix<double>& epsilon)const
 {
  // Get the Poisson ratio
  double nu(get_nu());
   
  // Now construct the Stress
  sigma(0,0) = (epsilon(0,0) + nu*epsilon(1,1)) / (1.0 - nu*nu);
  sigma(1,1) = (epsilon(1,1) + nu*epsilon(0,0)) / (1.0 - nu*nu);
  sigma(0,1) = epsilon(0,1) / (1.0 + nu);
  sigma(1,0) = sigma(0,1);
 }

 /// Get the principal stresses from the stress tensor
 void get_principal_stresses(const DenseMatrix<double>& sigma,
		Vector<double>& eigenvals,
		DenseMatrix<double>& eigenvecs)const
 {
  // Ensure that our eigenvectors are the right size
  eigenvals.resize(2);
  eigenvecs.resize(2);
  
  // Store the axial and shear stresses
  double s00 = sigma(0,0);
  double s01 = sigma(0,1);
  double s11 = sigma(1,1);

  // Calculate the principal stress magnitudes
  eigenvals[0] =
   0.5 * ( (s00 + s11) + sqrt((s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01)) );
  eigenvals[1] =
   0.5 * ( (s00 + s11) - sqrt((s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01)) );

  // Handle the shear free case
  if(s01==0.0)
   {
    eigenvecs(0,0)=1.0;
    eigenvecs(1,0)=0.0;
    eigenvecs(0,1)=0.0;
    eigenvecs(1,1)=1.0;
   }
  
  else
   {
    // TODO: better (more general) sign choice for streamlines

    // For max eval we choose y-positive evecs (suited to swelling sheet
    // problem)
    double sign = (eigenvals[0]-s00<0.0) ? -1.0 : 1.0;
    // Calculate the normalised principal stress direction for eigenvals[0]
    eigenvecs(0,0) =
     sign * (s01 / sqrt(s01*s01 + (eigenvals[0]-s00)*(eigenvals[0]-s00)));
    eigenvecs(1,0) =
     sign * ((eigenvals[0]-s00) / sqrt(s01*s01 + (eigenvals[0]-s00)*(eigenvals[0]-s00)));

    // For min eval we choose x-positive evecs (suited to swelling sheet
    // problem)
    sign = (s01<0.0) ? -1.0 : 1.0;
    // Calculate the normalised principal stress direction for eigenvals[1]
    eigenvecs(0,1) =
     sign * (s01 / sqrt(s01*s01 + (eigenvals[1]-s00)*(eigenvals[1]-s00)));
    eigenvecs(1,1) =
     sign * ((eigenvals[1]-s00) / sqrt(s01*s01 + (eigenvals[1]-s00)*(eigenvals[1]-s00)));
   }
 }

 /// Access function: Pointer to pressure function
 PressureFctPt& pressure_fct_pt() {return Pressure_fct_pt;}

 /// Access function: Pointer to pressure function. Const version
 PressureFctPt pressure_fct_pt() const {return Pressure_fct_pt;}

 /// Access function: Pointer to in plane forcing function
 InPlaneForcingFctPt& in_plane_forcing_fct_pt()
  {return In_plane_forcing_fct_pt;}

 /// Access function: Pointer to in plane forcing function. Const version
 InPlaneForcingFctPt in_plane_forcing_fct_pt() const
  {return In_plane_forcing_fct_pt;}

 /// Access function: Pointer to swelling function
 SwellingFctPt& swelling_fct_pt() {return Swelling_fct_pt;}

 /// Access function: Pointer to swelling function. Const version
 SwellingFctPt swelling_fct_pt() const {return Swelling_fct_pt;}

 /// Access function: Pointer to error metric function
 ErrorMetricFctPt& error_metric_fct_pt() {return Error_metric_fct_pt;}

 /// Access function: Pointer to multiple error metric function
 MultipleErrorMetricFctPt& multiple_error_metric_fct_pt() 
  {return Multiple_error_metric_fct_pt;}

 /// Access function: Pointer to error metric function function
 ErrorMetricFctPt error_metric_fct_pt() const {return Error_metric_fct_pt;}

 /// Access function: Pointer to multiple error metric function
 MultipleErrorMetricFctPt multiple_error_metric_fct_pt() const 
  {return Multiple_error_metric_fct_pt;}

 ///Access function to the Poisson ratio.
 const double*& nu_pt() {return Nu_pt;}

 ///Access function to the Poisson ratio (const version)
 const double& get_nu() const {return *Nu_pt;}

 ///Access function to the dampening coefficient.
 const double*& mu_pt() {return Mu_pt;}

 ///Access function to the dampening coefficient (const version)
 const double& get_mu() const {return *Mu_pt;}

 /// Get the kth dof type at internal point l at time t(=0)
 virtual double get_w_bubble_dof(const unsigned& l,
				 const unsigned& k,
				 const unsigned& t = 0) const =0;

 /// Get the kth equation at internal point l
 virtual int local_w_bubble_equation(const unsigned& l, const unsigned& k)
   const =0;

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline virtual void get_pressure_foeppl_von_karman(const unsigned& ipt,
                                        const Vector<double>& x,
                                        double& pressure) const
  {
   //If no pressure function has been set, return zero
   if(Pressure_fct_pt==0)
    {
     pressure = 0.0;
    }
   else
    {
     // Get pressure strength
     (*Pressure_fct_pt)(x,pressure);
    }
  }

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline virtual void get_in_plane_forcing_foeppl_von_karman(const unsigned& ipt,
                                        const Vector<double>& x,
                                        Vector<double>& pressure) const
  {
   //In plane is same as DIM of problem (2)
   pressure.resize(this->dim());
   //If no pressure function has been set, return zero
   if(In_plane_forcing_fct_pt==0)
    {
     pressure[0] = 0.0;
     pressure[1] = 0.0; 
    }
   else
    {
     // Get pressure strength
     (*In_plane_forcing_fct_pt)(x,pressure);
    }
  }

 /// Get swelling at (Eulerian) position x. This function is
 /// virtual to allow overloading.
 inline virtual void get_swelling_foeppl_von_karman(const unsigned& ipt,
                                        const Vector<double>& x,
                                        double& swelling) const
  {
   //If no swelling function has been set, return zero
   if(Swelling_fct_pt==0)
    {
     swelling = 0.0;
    }
   else
    {
     // Get swelling magnitude
     (*Swelling_fct_pt)(x,swelling);
    }
  }
 
 /// \short Calculate the elastic energy of the element and return it as a
 /// double.
 virtual Vector<double> element_elastic_and_kinetic_energy();

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_foeppl_von_karman(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }


 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_foeppl_von_karman(residuals,jacobian,1);
  }

 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {
   //Call fill in Jacobian 
   fill_in_contribution_to_jacobian(residuals,jacobian);
   // There is no mass matrix: we will just want J w = 0

   // -- COPIED FROM DISPLACMENT FVK EQUATIONS --
   // Dummy diagonal (won't result in global unit matrix but
   // doesn't matter for zero eigenvalue/eigenvector
   unsigned ndof=mass_matrix.nrow();
   for (unsigned i=0;i<ndof;i++)
    {
     mass_matrix(i,i)+=1.0;
    }
  }

 /// \short Return FE representation of dw/dt at local coordinate s
 virtual inline Vector<double> interpolated_dwdt_foeppl_von_karman(const Vector<double>& s) const
 {
  //Find out how many nodes positional dofs there are
  unsigned n_basis_type = nnodal_basis_type();
  // Find the internal dofs
  const unsigned n_b_position_type = nbubble_basis_type();
  //Find out how many nodes there are
  const unsigned n_w_node = nnode_outofplane();
  //Find out how many internal points there are
  const unsigned n_b_node = nbubble_basis();
  //Get the index at which the unknown is stored
  const unsigned w_nodal_index = w_nodal_index_foeppl_von_karman();
  //Get the index at which the unknown is stored
  const unsigned w_bubble_index = 1; //_foeppl_von_karman();
  // ^^^ THIS NEEDS TO BE FIXED, NOT GENERAL ENOUGH.
  // (1) should be changed to right index
  
  //Local c1-shape funtion
  Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
   psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
  
  //Calculate values of unknown
  Vector<double> interpolated_w(1,0.0);
  Vector<double> interpolated_dwdt(1,0.0);
    
  //Call the derivatives of the shape and test functions for the unknown
  shape_and_test_foeppl_von_karman(s,
				   psi, psi_b,
				   test, test_b);

  // Loop over nodes
  for(unsigned l=0;l<n_w_node;l++)
   {
    TimeStepper* timestepper_pt = this->node_pt(l)->time_stepper_pt();
    unsigned n_time = 0;
    if( !(timestepper_pt->is_steady()) )
     {
      n_time=timestepper_pt->ntstorage();
     }
    for(unsigned k=0;k<n_basis_type;k++)
     {
      // Add the contributions to the time derivative from node l, type k
      double dwdt_value = 0.0;
      for(unsigned t=0; t<n_time; t++)
       {
	dwdt_value +=
	 timestepper_pt->weight(1,t)
	 * this->raw_nodal_value(t, l, w_nodal_index+k);
       }
      interpolated_dwdt[0] += dwdt_value * psi(l,k);
     }
   }

  // Loop over internal dofs
  for(unsigned l=0;l<nbubble_basis();l++)
   {
    TimeStepper* timestepper_pt = internal_data_pt(w_bubble_index) // :(
     ->time_stepper_pt();
    unsigned n_time = 0;
    if( !(timestepper_pt->is_steady()) )
     {
      n_time=timestepper_pt->ntstorage();
     }
    for(unsigned k=0;k<n_b_position_type;k++)
     {
      // Add the contributions to the time derivative
      double dwdt_value = 0.0;
      for(unsigned t=0; t<n_time; t++)
       {
	dwdt_value +=
	 timestepper_pt->weight(1,t) * get_w_bubble_dof(l,k,t);
       }
      interpolated_dwdt[0] += dwdt_value * psi_b(l,k);
     }
   }
  // Finally return the (single valued) vector of interpolated velocity.
  return interpolated_dwdt;
 }
 
 /// \short Return FE representation of unknown values u(s)
 /// at local coordinate s
 virtual inline Vector<double> interpolated_u_foeppl_von_karman(const Vector<double> &s) const
  {
   //Find number of position dofs
   const unsigned n_basis_type = nnodal_basis_type();
   // Find the internal dofs
   const unsigned n_b_position_type = nbubble_basis_type();
   //Find out how many nodes there are
   const unsigned n_w_node = nnode_outofplane();
   //Find out how many nodes there are
   const unsigned n_u_node = nnode_inplane();
   //Find out how many internal points there are
   const unsigned n_b_node = nbubble_basis();
   //Get the index at which the unknown is stored
   const unsigned u_nodal_index = u_nodal_index_foeppl_von_karman();
   //Get the index at which the unknown is stored
   const unsigned w_nodal_index = w_nodal_index_foeppl_von_karman();
   //Get the index at which the unknown is stored
   const unsigned w_bubble_index = 0; //w_index_foeppl_von_karman();
   // Number of unique second deriv. will be the DIMth triangle number
   const unsigned n_second_deriv = (this->dim()*(this->dim()+1))/2;

   //Local c1-shape funtion
   Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
    psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
   
   DShape dpsi_dxi(n_w_node,n_basis_type,this->dim()),dtest_dxi(n_w_node,n_basis_type,this->dim()),
    dpsi_b_dxi(n_b_node,n_b_position_type,this->dim()),dtest_b_dxi(n_b_node,n_b_position_type,this->dim()),
    d2psi_dxi2(n_w_node,n_basis_type,n_second_deriv), d2test_dxi2(n_w_node,n_basis_type,n_second_deriv),
    d2psi_b_dxi2(n_b_node,n_b_position_type,n_second_deriv), d2test_b_dxi2(n_b_node,n_b_position_type,n_second_deriv);
   
   // In--plane dofs
   Shape psi_u(n_u_node);
   DShape dpsi_u(n_u_node,this->dim());
   Shape test_u(n_u_node);
   DShape dtest_u(n_u_node,this->dim());

   // Number of in-plane displacement fields is equal to dim + dim^2
   const unsigned n_u_fields = 6;// DIM;
   // Number of out-of-plane displacement fields is equal 1 (value of deflection_
   // plus first and second deriv
   const unsigned n_w_fields = 6;//1+DIM+n_second_deriv;
   //Initialise value of u
   Vector<double> interpolated_u(n_u_fields+n_w_fields,0.0);
   //Find values of c1-shape function
   d2shape_and_d2test_eulerian_foeppl_von_karman(s,psi,psi_b,dpsi_dxi,dpsi_b_dxi,
    d2psi_dxi2,d2psi_b_dxi2,test,test_b,dtest_dxi,dtest_b_dxi,d2test_dxi2,
    d2test_b_dxi2);

   // Get shape and test
   dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,psi_u,dpsi_u,test_u,dtest_u);
   //Interpolated unknown
   for(unsigned l=0;l<n_w_node;l++)
   {
    for(unsigned k=0;k<n_basis_type;k++)
     {
     // u_3
     interpolated_u[0] += this->nodal_value(l,w_nodal_index+k)*psi(l,k);
     // d_u_3_dx_alpha
     for(unsigned alpha=0;alpha<this->dim();++alpha)
      {interpolated_u[1+alpha] += this->nodal_value(l,w_nodal_index+k)*dpsi_dxi(l,k,alpha);}
     // d2_u_3_dx_alpha dx_beta
     for(unsigned alphabeta=0;alphabeta<n_second_deriv;++alphabeta)
      {
      interpolated_u[this->dim()+1+alphabeta] += this->nodal_value(l,w_nodal_index+k)
         *d2psi_dxi2(l,k,alphabeta);
      }
     }
   }

   // Bubble dofs
   for(unsigned l=0;l<nbubble_basis();l++)
    {
    for(unsigned k=0;k<nbubble_basis_type();k++)
     {
      double u_value = get_w_bubble_dof(l,w_bubble_index+k);
      // u_3
      interpolated_u[0] += u_value * psi_b(l,w_bubble_index+k);
      // d_u_3_dx_alpha
      for(unsigned alpha=0;alpha<this->dim();++alpha)
       { interpolated_u[1+alpha] += u_value*dpsi_b_dxi(l,w_bubble_index+k,alpha); }
      // d2_u_3_dx_alpha dx_beta
      for(unsigned alphabeta=0;alphabeta<n_second_deriv;++alphabeta)
       {
	interpolated_u[this->dim()+1+alphabeta] += u_value*d2psi_b_dxi2(l,w_bubble_index+k,alphabeta);
       }
     }
    }
   // Now for the displacement 
   for(unsigned l=0; l< n_u_node;++l)
    {
     // Now for the two in--plane displacements
     interpolated_u[6] += this->nodal_value(l,u_nodal_index+0)*psi_u(l);
     interpolated_u[7] += this->nodal_value(l,u_nodal_index+1)*psi_u(l);
     // Also output the in--plane displacement derivatives
     for(unsigned i=0; i<this->dim(); ++i)
      {
       interpolated_u[8 +i] += this->nodal_value(l,u_nodal_index+0)*dpsi_u(l,i);
       interpolated_u[10+i] += this->nodal_value(l,u_nodal_index+1)*dpsi_u(l,i);
      }
    }
   return(interpolated_u);
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();

protected:
 /// Pure virtual interface to the basis for the in--plane displacement 
 virtual void shape_u(const Vector<double> &s,  Shape &psi) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  DShape &d2psi_dx2,DShape& d2psi_b_dx2,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx,
  DShape &d2test_dx2,DShape& d2test_b_dx2) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx) const=0;

 /// \short Shape/test functions at local coordinate s
 virtual void shape_and_test_foeppl_von_karman(const Vector<double> &s,
  Shape &psi, Shape& psi_b, Shape &test, Shape& test_b) const=0;
 
 /// \short in--plane Shape/test functions at local coordinate s
 virtual double dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s,
  Shape &psi,DShape &dpsidx, Shape &test,DShape &dtestdx) const=0;

 /// \short Compute element residual Vector only (if flag=and/or element
 /// Jacobian matrix
 virtual void fill_in_generic_residual_contribution_foeppl_von_karman(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);

 /// Pointer to pressure function:
 PressureFctPt Pressure_fct_pt;

 /// Pointer to in plane forcing function (i.e. the shear force applied to
 /// the face)
 InPlaneForcingFctPt In_plane_forcing_fct_pt;
 
 /// Pointer to swelling function:
 SwellingFctPt Swelling_fct_pt; 

 /// Pointer to Poisson ratio, which this element cannot modify
 const double* Nu_pt;

 /// Pointer to the dampening coefficient, which this element cannot modify
 const double* Mu_pt;

 /// Pointer to global eta
 const double *Eta_pt;

 /// Default value for physical constant: Poisson ratio. 
 static const double Default_Nu_Value;

 /// Default value for constant: dampening coefficient. 
 static const double Default_Mu_Value;

 /// Default eta value so that we use 'natural' nondim and have no h dependence. 
 static const double Default_Eta_Value;

 /// Pointer to error metric
 ErrorMetricFctPt Error_metric_fct_pt;

 /// Pointer to error metric when we want multiple errors
 MultipleErrorMetricFctPt Multiple_error_metric_fct_pt;


 protected:

 /// Pointer to precomputed matrix that associates shape functions to monomials
 DenseMatrix<double> *Association_matrix_pt;
};


} //end namespace oomph
#endif

