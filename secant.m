## Copyright (C) 2022 caleb
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} secant (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: caleb <caleb@DESKTOP-SCRN63L>
## Created: 2022-04-18

function sol = secant (f,y1,y2,max_iters,tol)

y = zeros(max_iters,1);  %allocate memory for the solution guesses
%secant method requires 2 intial guesses: here, guess1 & guess2
y(1) = y1;  %initial guess "i" minus one (i-1)
y(2) = y2;  %initial guess "i"
i = 2;  %iteration counter; start @ 2 since secant requires y(i) - y(i-1)

%provide the initial residuals and convergences for the two initial guesses
res(1) = abs(f(y(1)));
res(2) = abs(f(y(2)));
conv(1) = 0;
conv(2) = abs(y(2) - y(1));

while i < max_iters
 
 %approximate the derivative of f(y)
  f_prime(i) = ( f(y(i)) - f(y(i-1)) ) ./ ( y(i) - y(i-1) );
  
  %find the next iteration's solution
  y(i+1) = y(i) - (f(y(i)) ./ f_prime(i));
  
  %calculate the residual and iterative convergence for this iteration
  res(i+1) = abs(f(y(i+1)));
  conv(i+1) = abs(y(i+1) - y(i));
  sol = y(i+1);

  %{
  if sol < 0
  printf("ERROR!! Solution reached a non-physical value at iteration %i.\n",i-1)
  break
  endif
  %}

  if (res(i+1) < tol) && (conv(i+1) < tol)
    break
  endif
  
  i+=1;
  
endwhile

printf("The solution is: %.10f meters.\n",sol)
printf("It took %i iterations to reach this solution using secant.\n",i-1)  
 
endfunction
