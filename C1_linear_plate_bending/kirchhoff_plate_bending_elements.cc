// Non--inline functions for the kirchhoff_plate_bending equations
#include "kirchhoff_plate_bending.h"

namespace oomph
{

//======================================================================
void  KirchhoffPlateBendingEquations::
fill_in_generic_residual_contribution_biharmonic(Vector<double> &residuals,
                                              DenseMatrix<double> &jacobian,
                                              const unsigned& flag)
{
 //Find out how many nodes there are
 const unsigned n_node = this->nnode();
 //Find out how many bubble nodes there are
 const unsigned n_b_node = nbubble_basis();
 //Find out how many nodes positional dofs there are
 unsigned n_basis_type = nnodal_basis_type();
 // Find the internal dofs
 const unsigned n_b_position_type = nbubble_basis_type();
 // Guaranteed to be an integer
 const unsigned n_2ndderiv = (dim()*(dim()+1))/2;

 //Get the Poisson ratio of the plate
 const double nu=get_nu();

 //Local c1-shape funtion
 Shape psi(n_node,n_basis_type),test(n_node,n_basis_type),
  psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
 
 DShape dpsi_dxi(n_node,n_basis_type,dim()),dtest_dxi(n_node,n_basis_type,dim()),
  dpsi_b_dxi(n_b_node,n_b_position_type,dim()),dtest_b_dxi(n_b_node,n_b_position_type,dim()),
  d2psi_dxi2(n_node,n_basis_type,n_2ndderiv), d2test_dxi2(n_node,n_basis_type,n_2ndderiv),
  d2psi_b_dxi2(n_b_node,n_b_position_type,n_2ndderiv), d2test_b_dxi2(n_b_node,n_b_position_type,n_2ndderiv);

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
   Vector<double> interpolated_u(1,0.0);
   DenseMatrix<double> interpolated_dudxi(1,dim(),0.0);
   DenseMatrix<double> interpolated_d2udxi2(1,n_2ndderiv,0.0);

   //Allocate and initialise to zero
   Vector<double> interpolated_x(dim(),0.0);
   Vector<double> s(dim());
   s[0] = this->integral_pt()->knot(ipt,0);
   s[1] = this->integral_pt()->knot(ipt,1);
   this->interpolated_x(s,interpolated_x);

   //Call the derivatives of the shape and test functions for the unknown
   double J = d2shape_and_d2test_eulerian_biharmonic(s,
    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate function value and derivatives
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++)
    {
     for(unsigned k=0;k<n_basis_type;k++)
      {
       //Get the nodal value of the unknown
       double u_value = this->raw_nodal_value(l,k);
       interpolated_u[0] += u_value*psi(l,k);
       // Loop over directions
       for(unsigned j=0;j<dim();j++)
        {
         interpolated_dudxi(0,j) += u_value*dpsi_dxi(l,k,j);
        }
       for(unsigned j=0;j<n_2ndderiv;j++)
        {
          interpolated_d2udxi2(0,j) += u_value*d2psi_dxi2(l,k,j);
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
      interpolated_u[0] += u_value*psi_b(l,k);
      // Loop over directions
      for(unsigned j=0;j<dim();j++)
       {
        interpolated_dudxi(0,j) += u_value*dpsi_b_dxi(l,k,j);
       }
      for(unsigned j=0;j<n_2ndderiv;j++)
       {
        interpolated_d2udxi2(0,j) += u_value*d2psi_b_dxi2(l,k,j);
       }
     }
    }

   //Get pressure function
   //-------------------
   double  pressure;
   get_pressure_biharmonic(ipt,interpolated_x,pressure);

   // Loop over the nodal test functions
   for(unsigned l=0;l<n_node;l++)
    {
     for(unsigned k=0;k<n_basis_type;k++)
      {
       //Get the local equation
       local_eqn = this->nodal_local_eqn(l,k);
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
           residuals[local_eqn] += (1-nu)*interpolated_d2udxi2(0,alpha+beta)
                                  *d2test_dxi2(l,k,alpha+beta)*W;
           // w_{,\alpha\alpha} \delta \kappa_\beta\beta
           residuals[local_eqn] += (nu)*interpolated_d2udxi2(0,beta+beta)
                                  *d2test_dxi2(l,k,alpha+alpha)*W;
           }
          }
         // Calculate the jacobian
         //-----------------------
         if(flag)
          {
           //Loop over the test functions again
           for(unsigned l2=0;l2<n_node;l2++)
            {
             // Loop over position dofs
             for(unsigned k2=0;k2<n_basis_type;k2++)
              {
               local_unknown = this->nodal_local_eqn(l2,k2);
               // If at a non-zero degree of freedom add in the entry
               if(local_unknown >= 0)
                {
                // Loop over dimensions 
                for(unsigned alpha=0;alpha<2;++alpha)
                 {
                 for(unsigned beta=0; beta<2;++beta)
                  {
                  // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                  jacobian(local_eqn,local_unknown) += (1-nu)
                     *d2psi_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
                  // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                  jacobian(local_eqn,local_unknown) += nu
                     *d2psi_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
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
                  // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                  jacobian(local_eqn,local_unknown) += (1-nu)
                     *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
                  // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                  jacobian(local_eqn,local_unknown) += nu
                     *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
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
           // w_{,\alpha\beta} \delta \kappa_\alpha\beta
           residuals[local_eqn] += (1-nu)*interpolated_d2udxi2(0,alpha+beta)
                                  *d2test_b_dxi2(l,k,alpha+beta)*W;
           // w_{,\alpha\alpha} \delta \kappa_\beta\beta
           residuals[local_eqn] += (nu)*interpolated_d2udxi2(0,beta+beta)
                                  *d2test_b_dxi2(l,k,alpha+alpha)*W;
           }
          }
         // Calculate the jacobian
         //-----------------------
         if(flag)
          {
           //Loop over the test functions again
           for(unsigned l2=0;l2<n_node;l2++)
            {
             // Loop over position dofs
             for(unsigned k2=0;k2<n_basis_type;k2++)
              {
               local_unknown = this->nodal_local_eqn(l2,k2);
               //If at a non-zero degree of freedom add in the entry
               if(local_unknown >= 0)
                {
                // Loop over dimensions 
                for(unsigned alpha=0;alpha<2;++alpha)
                 {
                 for(unsigned beta=0; beta<2;++beta)
                  {
                  // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                  jacobian(local_eqn,local_unknown) += (1-nu)
                     *d2psi_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
                  // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                  jacobian(local_eqn,local_unknown) += nu
                     *d2psi_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
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
                  // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                  jacobian(local_eqn,local_unknown) += (1-nu)
                     *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
                  // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                  jacobian(local_eqn,local_unknown) += nu
                     *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
                  }
                 }
                }
              }
            }
          } // End of flag
        }
      }
    }
  } // End of loop over integration points
}



//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
unsigned  KirchhoffPlateBendingEquations::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 // Return verdict
 if (passed)
  {
   return 0;
  }
 else
  {
   return 1;
  }

}

//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void  KirchhoffPlateBendingEquations::output(std::ostream &outfile,
                                    const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(dim()),x(dim());

 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);

 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 Vector<double> r(dim()+1);

 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   Vector<double> u(this->required_nvalue(0),0.0);
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   interpolated_u_biharmonic(s,u);

   // Get x position as Vector
   this->interpolated_x(s,x);
   for(unsigned i=0;i<dim();i++)
    {
     outfile << x[i] << " ";
    }

   // Loop for variables
   for(unsigned j=0;j<this->required_nvalue(0);j++)
    {
     outfile << u[j] << " " ;
    }

   outfile << std::endl;
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);
}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void  KirchhoffPlateBendingEquations::output(FILE* file_pt,
                                    const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(dim()), x(dim());;

 // Tecplot header info
 fprintf(file_pt,"%s",this->tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   Vector<double> u(this->required_nvalue(0),0.0);
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);

   // Get x position as Vector
   this->interpolated_x(s,x);

   for(unsigned i=0;i<dim();i++)
    {
     fprintf(file_pt,"%g ",x[i]);
    }
   interpolated_u_biharmonic(s,u);
   fprintf(file_pt,"%g \n",u[0]);//interpolated_u_poisson(s));
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(file_pt,nplot);
}



//======================================================================
 /// Output exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
void KirchhoffPlateBendingEquations::output_fct(std::ostream &outfile,
                                       const unsigned &nplot,
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(dim());

  // Vector for coordintes
  Vector<double> x(dim());

 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);

 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);

 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {

   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);

   // Get x position as Vector
   this->interpolated_x(s,x);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   //Output x,y,...,u_exact
   for(unsigned i=0;i<dim();i++)
    {
     outfile << x[i] << " ";
    }
   // Loop over variables
   for(unsigned j=0;j<this->required_nvalue(0);j++)
   {
    outfile << exact_soln[j] << " ";
   }
   outfile <<  std::endl;
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);
}




//======================================================================
 /// Validate against exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
//======================================================================
void KirchhoffPlateBendingEquations::compute_error(std::ostream &outfile,
                                          FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                                          double& error, double& norm)
{
 // Initialise
 error=0.0;
 norm=0.0;
 //Find out how many nodes there are
 const unsigned n_node = this->nnode();
 //Find out how many bubble nodes there are
 const unsigned n_b_node = nbubble_basis();
 //Find out how many nodes positional dofs there are
 unsigned n_basis_type = nnodal_basis_type();
 // Find the internal dofs
 const unsigned n_b_position_type = nbubble_basis_type();
 // Guaranteed to be an integer
 const unsigned n_2ndderiv = (dim()*(dim()+1))/2;
 //Local c1-shape funtion
 Shape psi(n_node,n_basis_type),test(n_node,n_basis_type),
  psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
 
 DShape dpsi_dxi(n_node,n_basis_type,dim()),dtest_dxi(n_node,n_basis_type,dim()),
  dpsi_b_dxi(n_b_node,n_b_position_type,dim()),dtest_b_dxi(n_b_node,n_b_position_type,dim()),
  d2psi_dxi2(n_node,n_basis_type,n_2ndderiv), d2test_dxi2(n_node,n_basis_type,n_2ndderiv),
  d2psi_b_dxi2(n_b_node,n_b_position_type,n_2ndderiv), d2test_b_dxi2(n_b_node,n_b_position_type,n_2ndderiv);

 //Vector of local coordinates
 Vector<double> s(dim());

 // Vector for coordintes
 Vector<double> x(dim());

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();

 // Tecplot
 //outfile << "ZONE" << std::endl;

 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<dim();i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J;
   // double Jlin = this->J_eulerian1(s);// Nope
   {
   J=this-> d2shape_and_d2test_eulerian_biharmonic(s,
    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);
   }
   //Premultiply the weights and the Jacobian
   double W = w*J;
   // double Wlin = w*Jlin;

   // Get x position as Vector
   this->interpolated_x(s,x);

   // Get FE function value
   Vector<double> u_fe(this->required_nvalue(0),0.0);
   interpolated_u_biharmonic(s,u_fe);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   //Output x,y,...,error
   for(unsigned i=0;i<dim();i++)
    {
     outfile << x[i] << " ";
    }
   for(unsigned ii=0;ii<this->required_nvalue(0);ii++)
    {
     outfile << exact_soln[ii] << " " << exact_soln[ii]-u_fe[ii] << " ";
    }
   outfile << std::endl;

   // Loop over variables
   double tmp1 = 0.0, tmp2 =0.0;
   for(unsigned ii=0;ii<1;ii++)
    {
     // Add to error and norm
     tmp1 = (exact_soln[ii]*exact_soln[ii]*W);
     tmp2 = ((exact_soln[ii]-u_fe[ii])*(exact_soln[ii]-u_fe[ii])*W);
     norm += tmp1;
     error += tmp2;
    }
  } //End of loop over integration pts
}

}
