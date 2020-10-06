# Option pricing using SSVI and local volatility model
A group project for Object Oriented Programming with Applications course at the University of Edinburgh taken for Computational Applied Mathematics MSc.
The C# program is focused on the local volatility model and the SSVI parametrization of the volatility surface. The program explores aspects of
option pricing using the SSVI and local volatility model. The program includes the following functions:

1. Call and Put prices from SSVI: Using the SSVI parametrisation class, the program obtains correct implied volatilities using Black–Scholes formula with
this implied volatility to price European put and call options.

2. Dupire SSVI Call/Put with Monte Carlo: Prices European put and call options in the SSVI local volatility framework using the Monte Carlo algorithm.

3. Pricing Asian arithmetic option: Prices Asian arithmetic call / put options using the Monte Carlo algorithm. 

4. Pricing lookback option: Prices a “lookback option” using the Monte Carlo algorithm.

# Contributors
Isabell Linde<br/>
Yanxin Fu<br/>
Yuexin Fui<br/>

# References
[1] Gawlikowicz, W. and Siska, D. A note on local volatility in option pricing<br/>
[2] Jim Gatheral & Antoine Jacquier (2014) Arbitrage-free SVI volatility surfaces,
Quantitative Finance, 14:1, 59-71, DOI: 10.1080/14697688.2013.819986<br/>
[3] Derman, Emanuel & Kani, Iraj. (1994). Riding on a Smile. Risk. 7. <br/>
[4] Gyöngy, I. Mimicking the one-dimensional marginal distributions of processes having an ito differential. Probab. Th. Rel. Fields 71, 501–516 (1986). https://doi.org/10.1007/BF00699039 <br/>
