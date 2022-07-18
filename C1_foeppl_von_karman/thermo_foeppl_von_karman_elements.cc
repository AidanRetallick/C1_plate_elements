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
// Non--inline functions for the foeppl_von_karman equations
#include "thermo_foeppl_von_karman_elements.h"

namespace oomph
{
//======================================================================
/// Fill in generic residual and jacobian contribution 
//======================================================================
void  ThermoFoepplVonKarmanEquations::
fill_in_generic_residual_contribution_thermo_foeppl_von_karman(
 Vector<double> &residuals,
 DenseMatrix<double> &jacobian,
 const unsigned& flag)
{
 //Find out how many nodes there are
 const unsigned n_u_node = nnode_inplane();
 const unsigned n_w_node = nnode_outofplane(); 
 //Find out how many bubble nodes there are
 const unsigned n_b_node = nbubble_basis();
 //Find out how many nodes positional dofs there are
 unsigned n_basis_type = nnodal_basis_type();
 // Find the internal dofs
 const unsigned n_b_position_type = nbubble_basis_type();

 //Get the Poisson ratio of the plate
 const double nu=get_nu();

 //Set up memory for the shape and test functions
 Shape psi_u(n_u_node), test_u(n_u_node);
 DShape dpsi_udxi(n_u_node,this->dim()), dtest_udxi(n_u_node,this->dim());

 //Local c1-shape funtion
 Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
  psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
 
 DShape dpsi_dxi(n_w_node,n_basis_type,this->dim()),dtest_dxi(n_w_node,n_basis_type,this->dim()),
  dpsi_b_dxi(n_b_node,n_b_position_type,this->dim()),dtest_b_dxi(n_b_node,n_b_position_type,this->dim()),
  d2psi_dxi2(n_w_node,n_basis_type,3), d2test_dxi2(n_w_node,n_basis_type,3),
  d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

 //Set the value of n_intpt
 const unsigned n_intpt = this->integral_pt()->nweight();
 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   //Calculate values of unknown
   Vector<double> interpolated_w(1,0.0);
   DenseMatrix<double> interpolated_dwdxi(1,this->dim(),0.0);
   DenseMatrix<double> interpolated_d2wdxi2(1,3,0.0);
   
   //Calculate values of unknown
   Vector<double> interpolated_u(2,0.0);
   DenseMatrix<double> interpolated_duidxj(2,this->dim(),0.0);

   //Allocate and initialise to zero
   Vector<double> interp_x(this->dim(),0.0);
   Vector<double> s(this->dim());
   s[0] = this->integral_pt()->knot(ipt,0);
   s[1] = this->integral_pt()->knot(ipt,1);
   interpolated_x(s,interp_x);
   //Call the derivatives of the shape and test functions for the unknown
   double J = d2shape_and_d2test_eulerian_foeppl_von_karman(s,
    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);

   dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,psi_u,dpsi_udxi,test_u,dtest_udxi);
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate function value and derivatives
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_w_node;l++)
    {
     for(unsigned k=0;k<n_basis_type;k++)
      {
       //Get the nodal value of the unknown
       double w_value = this->raw_nodal_value(l,k+2);
       interpolated_w[0] += w_value*psi(l,k);
       // Loop over directions
       for(unsigned j=0;j<this->dim();j++)
        {
         interpolated_dwdxi(0,j) += w_value*dpsi_dxi(l,k,j);
        }
       for(unsigned j=0;j<3;j++)
        {
          interpolated_d2wdxi2(0,j) += w_value*d2psi_dxi2(l,k,j);
        }
      }
    }

   // Loop over internal dofs
   for(unsigned l=0;l<nbubble_basis();l++)
   {
    for(unsigned k=0;k<n_b_position_type;k++)
     {
      //Get the nodal value of the unknown
      double u_value = get_w_bubble_dof(l,k);
      interpolated_w[0] += u_value*psi_b(l,k);
      // Loop over directions
      for(unsigned j=0;j<this->dim();j++)
       {
        interpolated_dwdxi(0,j) += u_value*dpsi_b_dxi(l,k,j);
       }
      for(unsigned j=0;j<3;j++)
       {
        interpolated_d2wdxi2(0,j) += u_value*d2psi_b_dxi2(l,k,j);
       }
     }
    }

   // Loop over nodes
   for(unsigned l=0;l<n_u_node;l++)
    {
     for(unsigned alpha=0;alpha<2;alpha++)
      {
       double u_value = this->raw_nodal_value(l,alpha);
       interpolated_u[alpha] += u_value*psi_u(l);
  
        // Loop over directions
        for(unsigned beta=0; beta<this->dim(); beta++)
         {
          interpolated_duidxj(alpha,beta) += u_value*dpsi_udxi(l,beta);
         }
      }
     }
   
   //Get the temperature function
   //--------------------
   double c_swell(0.0);
   get_swelling_foeppl_von_karman(ipt,interp_x,c_swell);

   // Get the stress
   DenseMatrix<double> sigma(2,2,0.0);
   get_sigma(sigma,interpolated_duidxj, interpolated_dwdxi, c_swell); 
   
    // Check 
//    Vector<double> interpolated_w(6,0.0)pli:
//    interpolated_w=this-> interpolated_w_foeppl_von_karman(s);
//   interpolated_w[0]=interpolated_w[0];
//   interpolated_dwdxi(0,0)=interpolated_w[1];
//   interpolated_dwdxi(0,1)=interpolated_w[2];
//   interpolated_d2wdxi2(0,0)=interpolated_w[3];
//   interpolated_d2wdxi2(0,1)=interpolated_w[4];
//   interpolated_d2wdxi2(0,2)=interpolated_w[5];
   //Get pressure function
   //-------------------
   double pressure(0.0);
   Vector<double> pressure_gradient(2,0.0);
   get_pressure_foeppl_von_karman(ipt,interp_x,pressure);
   get_in_plane_forcing_foeppl_von_karman(ipt,interp_x,pressure_gradient);

   // Loop over the nodal test functions
   for(unsigned l=0;l<n_w_node;l++)
    {
    for(unsigned k=0;k<n_basis_type;k++)
     {
     //Get the local equation
     local_eqn = this->nodal_local_eqn(l,k+2);
     //IF it's not a boundary condition
     if(local_eqn >= 0)
      {
      // Add body force/pressure term here
      residuals[local_eqn] -= pressure*test(l,k)*W;
      for(unsigned alpha=0;alpha<2;++alpha)
       {
       for(unsigned beta=0; beta<2;++beta)
        {
        // w_{,\alpha\beta} \delta \kappa_\alpha\beta
        residuals[local_eqn] += (1-nu)
         *interpolated_d2wdxi2(0,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
        // w_{,\alpha\alpha} \delta \kappa_\beta\beta
        residuals[local_eqn] += (nu)
         *interpolated_d2wdxi2(0,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
        // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
        residuals[local_eqn] += eta()*sigma(alpha,beta)
         *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
        }
       }
      // Calculate the jacobian
      //-----------------------
      if(flag)
      {
       //Loop over the in--plane unknowns again
       for(unsigned ll=0;ll<n_u_node;ll++)
        {
        // Loop over displacements
        for(unsigned gamma=0;gamma<2;++gamma)
         {
         local_unknown = this->nodal_local_eqn(ll,gamma);
         // If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over displacements
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           // Loop over dimensions 
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Nonlinear Terms 
            if(alpha==beta)
             {
             jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)* 
                   interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(alpha==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
                   *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(beta==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
                   *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            }
           }
          }
         }
        }

       //Loop over the test functions again
       for(unsigned l2=0;l2<n_w_node;l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_basis_type;k2++)
         {
          local_unknown = this->nodal_local_eqn(l2,k2+2);
          // If at a non-zero degree of freedom add in the entry
          if(local_unknown >= 0)
           {
           // Loop over dimensions 
           for(unsigned alpha=0;alpha<2;++alpha)
            {
            for(unsigned beta=0; beta<2;++beta)
             {
             // Linear Terms
             // w_{,\alpha\beta} \delta \kappa_\alpha\beta
             jacobian(local_eqn,local_unknown) += (1-nu)
                *d2psi_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
             // w_{,\alpha\alpha} \delta \kappa_\beta\beta
             jacobian(local_eqn,local_unknown) += nu
                *d2psi_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                        *dpsi_dxi(l2,k2,alpha)*dtest_dxi(l,k,beta)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,beta)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,beta)
               *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,beta)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
             }
            }
           }
         }
        }
       //Loop over the internal test functions
       for(unsigned l2=0;l2<nbubble_basis();l2++)
        {
        for(unsigned k2=0;k2<n_b_position_type;k2++)
         {
         local_unknown = local_w_bubble_equation(l2,k2);
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over dimensions 
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Linear Terms 
            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
            jacobian(local_eqn,local_unknown) += (1-nu)
               *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
            // w_{,\alpha\alpha} \delta \kappa_\beta\beta
            jacobian(local_eqn,local_unknown) += nu
               *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_dxi(l,k,beta)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,beta)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
              *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
            }
           }
          }
         }
        }
       } // End of flag
      }
     }
    }

    // Loop over the internal test functions
    for(unsigned l=0;l<nbubble_basis();l++)
    {
    for(unsigned k=0;k<n_b_position_type;k++)
     {
     //Get the local equation
     local_eqn = local_w_bubble_equation(l,k);
     //IF it's not a boundary condition
     if(local_eqn >= 0)
      {
      // Add body force/pressure term here
      residuals[local_eqn] -= pressure*test_b(l,k)*W;

      for(unsigned alpha=0;alpha<2;++alpha)
       {
       for(unsigned beta=0; beta<2;++beta)
        {
        // The Linear Terms 
        // w_{,\alpha\beta} \delta \kappa_\alpha\beta
        residuals[local_eqn] += (1-nu)
          *interpolated_d2wdxi2(0,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
        // w_{,\alpha\alpha} \delta \kappa_\beta\beta
        residuals[local_eqn] += (nu)
          *interpolated_d2wdxi2(0,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
        // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
        residuals[local_eqn] += eta()*sigma(alpha,beta)
          *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
        }
       }
      // Calculate the jacobian
      //-----------------------
      if(flag)
       {
       //Loop over the in--plane unknowns again
       for(unsigned ll=0;ll<n_u_node;ll++)
        {
        // Loop over displacements
        for(unsigned gamma=0;gamma<2;++gamma)
         {
         local_unknown = this->nodal_local_eqn(ll,gamma);
         // If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over displacements
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           // Loop over dimensions 
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Nonlinear Terms 
            if(alpha==beta)
             {
             jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)*
                   interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(alpha==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
                   *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(beta==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
                   *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            }
           }
          }
         }

        }
       //Loop over the test functions again
       for(unsigned l2=0;l2<n_w_node;l2++)
        {
         // Loop over position dofs
         for(unsigned k2=0;k2<n_basis_type;k2++)
          {
          local_unknown = this->nodal_local_eqn(l2,k2+2);
          //If at a non-zero degree of freedom add in the entry
          if(local_unknown >= 0)
           {
           // Loop over dimensions 
           for(unsigned alpha=0;alpha<2;++alpha)
            {
            for(unsigned beta=0; beta<2;++beta)
             {
             // The Linear Terms 
             // w_{,\alpha\beta} \delta \kappa_\alpha\beta
             jacobian(local_eqn,local_unknown) += (1-nu)
                *d2psi_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
             // w_{,\alpha\alpha} \delta \kappa_\beta\beta
             jacobian(local_eqn,local_unknown) += nu
                *d2psi_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                        *dpsi_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,beta)
               *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,beta)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
             }
            }
           }
          }
         }
       //Loop over the test functions again
       for(unsigned l2=0;l2<nbubble_basis();l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_b_position_type;k2++)
         {
         local_unknown = local_w_bubble_equation(l2,k2);
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over dimensions 
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Linear Terms 
            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
            jacobian(local_eqn,local_unknown) += (1-nu)
               *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
            // w_{,\alpha\alpha} \delta \kappa_\beta\beta
            jacobian(local_eqn,local_unknown) += nu
               *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
              *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
            }
           }
          }
         }
        }
       } // End of flag
      }
     }
    } //end loop over nodes
  for(unsigned l=0;l<n_u_node;++l)
   {
    // Now loop over displacement equations
    for(unsigned alpha=0; alpha<2 ; ++alpha)
     {
     //Get the local equation
     local_eqn = this->nodal_local_eqn(l,alpha);
     //IF it's not a boundary condition
     if(local_eqn >= 0)
      {
      // Now add on the stress terms
      for(unsigned beta=0; beta<2;++beta)
        {
         // Variations in u_alpha
         residuals[local_eqn] += sigma(alpha,beta)*dtest_udxi(l,beta)*W;
        }
      // Add the forcing term
      residuals[local_eqn] += pressure_gradient[alpha]*psi_u(l)*W;

      // Now loop over Jacobian
      if(flag)
       {
       for(unsigned ll=0;ll<n_u_node;++ll)
        {
        // Now loop over displacement equations
        for(unsigned gamma=0; gamma<2 ; ++gamma)
         {
         //Get the local equation
         local_unknown = this->nodal_local_eqn(ll,gamma);
         // If the unknown is a degree of freedom 
         if(local_unknown>=0)
          {
          // Now add on the stress terms
          for(unsigned beta=0; beta<2;++beta)
           {
           // The Linear Terms 
           if(alpha==beta)
            {
            jacobian(local_eqn,local_unknown) += nu*dpsi_udxi(ll,gamma)
                  *dtest_udxi(l,beta)/(1-nu*nu)*W;
            }
           if(alpha==gamma)
            {
            jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,beta)
                  *dtest_udxi(l,beta)/(1-nu*nu)*W;
            }
           if(beta==gamma)
            {
            jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,alpha)
                  *dtest_udxi(l,beta)/(1-nu*nu)*W;
            }
           }// End loop beta
          }// End if local unknown
         }// End loop gamma dof
        }// End loop nodal dofs
       
       //Loop over the w nodal dofs
       for(unsigned l2=0;l2<n_w_node;l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_basis_type;k2++)
         {
         //If at a non-zero degree of freedom add in the entry
         //Get the local equation
         local_unknown = this->nodal_local_eqn(l2,k2+2);
         // If the unknown is a degree of freedom 
         if(local_unknown>=0)
          {
          // Now add on the stress terms
          for(unsigned beta=0; beta<2;++beta)
           {
           // The NonLinear Terms 
           jacobian(local_eqn,local_unknown) += nu*dpsi_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_dxi(l2,k2,alpha)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           }// End loop beta
          }// End if local unknown
         }// End loop position type dof
        }// End loop w nodal dofs
       
       //Loop over the w internal dofs
       for(unsigned l2=0;l2<nbubble_basis();l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_b_position_type;k2++)
         {
         //If at a non-zero degree of freedom add in the entry
         //Get the local equation
         local_unknown = local_w_bubble_equation(l2,k2);
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Now add on the stress terms
          for(unsigned beta=0; beta<2;++beta)
           {
           // The NonLinear Terms 
           jacobian(local_eqn,local_unknown) += nu*dpsi_b_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_b_dxi(l2,k2,alpha)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_b_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           } // End loop beta
          } // End if local_unknown
         } // End loop position type
        } // End loop internal dofs
       } // End if flag
      } // End if local eqn
     } // End loop over alpha
    } // End loop over nodes
  } // End of loop over integration points
}

}

