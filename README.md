# Notes
* All of the code herein has been developed by David Robinson and Aidan Retallick (University of 
  Manchester) from existing code within oomph-lib.
# sturdy-palm-tree
* Generalised Foppl von Karman Equations using C1 triangular elements packaged
  as a user_src/ 

### Quick command for find and replace option on all .cc and .h files
```
find ../ -regex ".*\.[c,h]c?" -exec sed -i 's/matchstr/replstr/g' {} \;
```
