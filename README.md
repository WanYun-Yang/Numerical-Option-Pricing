## Numerical Option Pricing

Ting-Wei Su, Chien-Hao Wu, Wan-Yun Yang <br>
Department of Mathematics, National Chiao Tung University <br>
1001 University Road, Hsinchu, Taiwan 300, ROC <br>
Jan 8, 2018 <br>

&emsp; In this project, we focused on the numerical performance of different options and used three prevailing discretization schemes to make an extensive analysis of the options considered. We want to solve the Black-Scholes equation in the following three ways to compare the performance and accuracy: forward Euler method, Crank-Nicolson method, and solve a parabolic equation transformed from the Black-Scholes equation. We regarded the solution obtained by a binomial tree as the exact solution.

### An Introduction to Black-Scholes Equation and Options
&emsp; This is the page of<a href="https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_equation" title="Title">
Black-Scholes equation</a> on Wikipedia. 
We are keen to examine four different options: <a href="https://www.investopedia.com/terms/e/europeanoption.asp" title="Title">
European call</a>, <a href="https://www.investopedia.com/terms/e/europeanoption.asp" title="Title">
European put</a>, <a href="https://www.investopedia.com/terms/a/americanoption.asp" title="Title">
American call</a>, <a href="https://www.investopedia.com/terms/a/americanoption.asp" title="Title">
American call</a>. 

### Three Ways to Solve the Black-Scholes Equation
* Forward Euler method(FD_OptionPricing.m). <br>
&emsp; A simple and explicit <a href="https://en.wikipedia.org/wiki/Finite_difference" title="Title">
 finite difference method</a> to numerically solve partial differential equations. Notice that it must satisfy <a href="https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition" title="Title">
a certain CFL condition</a>. This means it's conditionally stable. <br> <br>
* Direct Crank-Nicolson method(CN_OptionPricing.m). <br>
&emsp; An implicit finite difference method to numerically solve partial differential equations. It's unconditionally stable, but we have to solve a tridiagonal system each step. To speed up the calculation, we use <a href="https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm" title="Title">the Thomas algorithm</a> to deal with the tridiagonal system. <br> <br>
* Crank-Nicolson method for the equaiton under heat-equation transformation(HT_OptionPricing.m). <br>
&emsp; See <a href="https://fenix.tecnico.ulisboa.pt/downloadFile/395139424085/Extended%20Abstract.pdf" title="Title">Finite Differences Schemes for Pricing of European and American Options</a>, or check Mathematics of Financial Derivatives a Student Introduction, Wilmott, 1995, page 76-81.

### Error Analysis
&emsp; The error comes from two places:
* The upper bound for the spatial boundary <br>
&emsp; When we do numerical calculations, we need to set the maximum stock price S<sub>max</sub> and obtain the numerical solution in (0, T) X (0, S<sub>max</sub>). The solution converges as S<sub>max</sub> goes to infinity.
* Local truncation error (LTE) <br>
&emsp; The error occurs from truncating the Taylor expansion into finite sums. The solution converges as temporal and spatial steps go to zero.

&emsp; Combined the two reasons mentioned above, the maximum absolute error (take L-infinity-norm on the error vector) is about 10<sup>-6</sup>. Though we cannot obtain as accurate as possible, this result is quite adequate for practical purposes.

### Conclusion
&emsp;We use three ways to calculate option price and all of them converge. 
