## Numerical Option Pricing

<a href="https://www.linkedin.com/in/ting-wei-su-2b508614b/" title="Title">Ting-Wei Su</a>, <a href="https://www.linkedin.com/in/chien-hao-wu-252a77137/" title="Title">Chien-Hao Wu</a>, <a href="https://www.linkedin.com/in/wan-yun-yan-b627b2121/" title="Title">Wan-Yun Yang</a> <br>
<a href="http://www.math.nctu.edu.tw/event/" title="Title">Department of Mathematics</a>, <a href="http://www.nctu.edu.tw/" title="Title">National Chiao Tung University</a> <br>
1001 University Road, Hsinchu, Taiwan 300, ROC <br>
Jan 8, 2018 <br>

&emsp; In this project, we focused on the numerical performance of different options and used three common discretization schemes to make an extensive analysis of the options considered. We want to solve the Black-Scholes equation in the following three ways to compare the performance and accuracy: <a href="https://en.wikipedia.org/wiki/Euler_method" title="Title">forward Euler method</a>, <a href="https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method" title="Title">Crank-Nicolson method</a>, and solve a parabolic equation transformed from the Black-Scholes equation. We regarded the solution obtained by a binomial tree as the exact solution. For one of the schemes, the forward Euler method, we presented two main propositions concerning stability and convergence of the Black-Scholes equation. We also applied the Crank-Nicolson method to the original Black-Scholes equation and the one under heat-equation transformation. Finally, in order to correctly solve the equation, we gave detailed explanations for the practical modifications of the boundary conditions in both European and American options. The related materials for the project can be found in <a href="https://github.com/WanYun-Yang/Numerical-Option-Pricing/blob/master/Handout.pdf" title="Title">Handout.pdf</a>.

### An Introduction to Black-Scholes Equation and Options
&emsp; This is the page of<a href="https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_equation" title="Title">
Black-Scholes equation</a> on Wikipedia. 
We are keen to examine four different options: <a href="https://www.investopedia.com/terms/e/europeanoption.asp" title="Title">
European call</a>, <a href="https://www.investopedia.com/terms/e/europeanoption.asp" title="Title">
European put</a>, <a href="https://www.investopedia.com/terms/a/americanoption.asp" title="Title">
American call</a> and <a href="https://www.investopedia.com/terms/a/americanoption.asp" title="Title">
American call</a>. 

### Three Ways to Solve the Black-Scholes Equation
* Forward Euler method (<a href="https://github.com/WanYun-Yang/Numerical-Option-Pricing/blob/master/FD_OptionPricing.m" title="Title">FD_OptionPricing.m</a>) <br> 
&emsp; A simple and explicit finite difference method to numerically solve partial differential equations. Notice that it must satisfy <a href="https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition" title="Title">
a certain CFL condition</a>. This means it's conditionally stable. <br> <br>
* Direct Crank-Nicolson method (<a href="https://github.com/WanYun-Yang/Numerical-Option-Pricing/blob/master/CN_OptionPricing.m" title="Title">CN_OptionPricing.m</a>) <br>
&emsp; An implicit finite difference method to numerically solve partial differential equations. It's unconditionally stable, but we have to solve a tridiagonal system at each step in time marching. To speed up the calculation, we use the<a href="https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm" title="Title"> Thomas algorithm</a> to deal with the tridiagonal system. <br> <br>
* Crank-Nicolson method for the equaiton under heat-equation transformation (<a href="https://github.com/WanYun-Yang/Numerical-Option-Pricing/blob/master/HT_OptionPricing.m" title="Title">HT_OptionPricing.m</a>) <br>
&emsp; See <a href="https://fenix.tecnico.ulisboa.pt/downloadFile/395139424085/Extended%20Abstract.pdf" title="Title">Finite Differences Schemes for Pricing of European and American Options</a>, or check <a href="https://books.google.com.tw/books?id=VYVhnC3fIVEC&printsec=frontcover&dq=Mathematics+of+Financial+Derivatives+a+Student+Introduction&hl=zh-TW&sa=X&ved=0ahUKEwiXuqKWxu3YAhXGlJQKHfz1BgYQ6AEIJjAA#v=onepage&q=Mathematics%20of%20Financial%20Derivatives%20a%20Student%20Introduction&f=false" title="Title">The Mathematics of Financial Derivatives: A Student Introduction</a>, <a href="https://en.wikipedia.org/wiki/Paul_Wilmott" title="Title">Wilmott</a>, 1995, page 76-81.



### Error Analysis
&emsp; The error comes from two places:
* The upper bound for the spatial boundary <br>
&emsp; When we do numerical calculations, we need to set the maximum stock price S<sub>max</sub> and obtain the numerical solution in (0, T) x (0, S<sub>max</sub>). The solution converges as S<sub>max</sub> goes to infinity.
* Local truncation error (LTE) <br>
&emsp; The error occurs from truncating the <a href="https://en.wikipedia.org/wiki/Taylor_series" title="Title">Taylor expansion</a> into finite sums. The solution converges as temporal and spatial steps go to zero.

&emsp; Combined the two reasons mentioned above, the maximum absolute error (take <a href="http://mathworld.wolfram.com/L-Infinity-Norm.html" title="Title">L^infty-norm</a> on the error vector) is about 10<sup>-6</sup>. Though we cannot obtain the numerical solution as accurate as possible, this result is quite adequate for practical purposes.

### Main Results
&emsp; We present a performance ranking and a convergence analysis for each of the three schemes.
