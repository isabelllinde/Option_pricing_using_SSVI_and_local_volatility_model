using System;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace Project2
{
    class Gaussian
    {
        //Gaussian function for approximation for cumulative normal distribution
        static public double psi(double x)
        //The probability density function
        {
            double A = 1.0 / Math.Sqrt(2.0 * Math.PI);
            return A * Math.Exp(-x * x * 0.5);
        }

        // Approximation for cumulative normal distribution used in B-S Formula. Values from Page 34 of C# for Financial Markets Book.
        static public double N(double x)
        {
            System.Diagnostics.Debug.Assert(!Double.IsNaN(x) && !Double.IsInfinity(x),
                "Gaussian.N: x has to be a valid number.");
            double a1 = 0.4361836;
            double a2 = -0.1201676;
            double a3 = 0.9372980;
            double k = 1.0 / (1.0 + (0.33267 * x));
            if (x >= 0.0)
            {
                return 1.0 - psi(x) * (a1 * k + (a2 * k * k) + (a3 * k * k * k));
            }
            else
            {
                return 1.0 - N(-x);
            }
        }
    }

    public class SSVI
    {
        // Initializing the parameters
        // These are private since they are only used in methods within class SSVI
        private double alpha = Math.Sqrt(0.1);
        private double beta = 1;
        private double gamma = 0.7;
        private double eta = 0.2;
        private double rho = 0.3;

        // S and r are public because S and r used in other classe
        public double S = 100;
        public double r = 0.025;

        
        public double Theta(double T)
        {
            return Math.Pow(alpha, 2) * (Math.Exp(Math.Pow(beta, 2) * T) - 1); //Theta(t) in SSVI parametrization
        }

        public double Logmoney(double S, double K, double T)
        {
            return Math.Log(K * Math.Exp(-r * T) / S); //Log-moneyness k in SSVI parametrisation
        }

        public double Omega(double S, double K, double T) // Parametrization function defined with Logmoneyness
        {

            double Phi = eta / (Math.Pow(Theta(T), gamma) * Math.Pow(1 + Theta(T), 1 - gamma)); //Phi(theta) in SSVI parametrization

            double omega_ssvi = ((Theta(T) / 2) * (1 + rho * Phi * Logmoney(S, K, T) + Math.Sqrt(Math.Pow((Phi * Logmoney(S,K,T) + rho), 2) + Math.Pow((1 - rho), 2))));

            return omega_ssvi;

        }

        public double Omega2(double K, double T) // Parametrization function defined without Logmoneyness
        {

            double Phi = eta / (Math.Pow(Theta(T), gamma) * Math.Pow(1 + Theta(T), 1 - gamma)); //Phi(theta) in SSVI parametrisation

            double omega_ssvi = ((Theta(T) / 2) * (1 + rho * Phi + Math.Sqrt(Math.Pow((Phi + rho), 2) + Math.Pow((1 - rho), 2))));

            return omega_ssvi;

        }
        public double ImpliedVolatility_atmoney(double T) // Implied Volatility sigma from SSVI when option prices at the money (k=0)
        {
            double sigma = Math.Sqrt(Theta(T) / T);

            return sigma;
        }
        public double ImpliedVolatility(double S, double K, double T) // Implied Volatility sigma from SSVI
        {
            double sigma = Math.Sqrt(Omega(S, K, T) / T);

            return sigma;
        }
    }

    class BlackScholesFormula : SSVI //Create a class to calculate the Black-Schole forumla, inheriting parameters and methods from SSVI class
    {
        public double CalculateCallOptionPrice(double S, double K, double T) //Black Scholes formula for calculating call option price
        {
            double sigma = ImpliedVolatility(S, K, T);

            if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

            double d1 = (Math.Log(S / K) + (r + (sigma * sigma) / 2) * (T)) / (sigma * Math.Sqrt(T));
            double d2 = d1 - sigma * Math.Sqrt(T);
            return S * Gaussian.N(d1) - K * Math.Exp(-r * (T)) * Gaussian.N(d2);
        }

        public static double GetCallFromPutPrice(double S, double K, double r, double T, double putPrice)
        {
            return putPrice - Math.Exp(-r * T) * K + S;
        }

        public static double GetPutFromCallPrice(double S, double K, double r, double T, double callPrice)
        {
            return callPrice - S + K * Math.Exp(-r * T);
        }

        public double CalculatePutOptionPrice(double S, double K, double T) //Black Scholes formula for calculating put option price
        {
            double sigma = ImpliedVolatility(S, K, T);

            if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

            return GetPutFromCallPrice(S, K, r, T, CalculateCallOptionPrice(S, K, T));
        }


        //at money using different omega formula
        public double CalculateCallOptionPrice_atmoney(double S, double K, double T)
        {
            double sigma = ImpliedVolatility_atmoney(T);

            if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

            double d1 = (Math.Log(S / K) + (r + (sigma * sigma) / 2) * (T)) / (sigma * Math.Sqrt(T));
            double d2 = d1 - sigma * Math.Sqrt(T);
            return S * Gaussian.N(d1) - K * Math.Exp(-r * (T)) * Gaussian.N(d2);
        }

        public static double GetCallFromPutPrice_atmoney(double S, double K, double r, double T, double putPrice)
        {
            return putPrice - Math.Exp(-r * T) * K + S;
        }

        public static double GetPutFromCallPrice_atmoney(double S, double K, double r, double T, double callPrice)
        {
            return callPrice - S + K * Math.Exp(-r * T);
        }

        public double CalculatePutOptionPrice_atmoney(double S, double K, double T)
        {
            double sigma = ImpliedVolatility_atmoney(T);

            if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

            return GetPutFromCallPrice(S, K, r, T, CalculateCallOptionPrice_atmoney(S, K, T));
        }
    }

    class BlackScholesMonteCarlo : SSVI
    {
        public double SigmaDupSquared(double S, double K, double T)
        {
            if (S <= 0 || K <= 0 || T <= 0)
                throw new System.ArgumentException("Need S >0, K > 0 and T > 0.");

            double eps = 1e-6;
            double h = 1e-6;

            //Using central finite difference method to calculate the partial derivatives
            double partialk = (Omega2(Logmoney(S, K, T) + h, T) - Omega2((Logmoney(S, K, T) - h), T) / (2 * h)); //Approximating of partial derivative with respect to k at point (T,K)
            if (partialk > eps)
            {
                return (Omega(S, K, T + h) - Omega(S, K, T - h)) / (2 * h); //Approximating of the partial derivative with respect to t at the point (T, K)
            }
            else
            {
                return (Omega(S, K, T + h) - Omega(S, K, T - h)) / (2 * h * G(S, K, T)); //Approximation of the partial derivative with respect to t at the point (T, K) divided by function G
            }
        }

        public double G(double S, double K, double T) // function g(t,k) in the Dupire formula
        {
            // exception handling
            if (T <= 0 || K <= 0 || S <= 0)
                throw new ArgumentException("Need S > 0, K > 0 and T > 0.");

            double partialk = SigmaDupSquared(S, K, T);

            double h = 1e-6;
            double partialder1 = (Omega2(Logmoney(S, K, T) + h, T) - Omega2(Logmoney(S, K, T) - h, T) / (2 * h)); //Approximation of the first derivative of with respect to k at point (T,K)
            double partialder2 = (Omega2(Logmoney(S, K, T) + h, T) - 2 * Omega2(Logmoney(S, K, T), T) + Omega2(Logmoney(S, K, T) - h, T) / (h * h)); //Approximation of the second derivative with respect to k at point (T,K)
            return Math.Pow(1 - (Logmoney(S, K, T) * partialk) / (2 * Omega(S, K, T)), 2) - 0.25 * partialder1 * partialder2 * (0.25 + 1 / Omega(S, K, T)) + 0.5 * partialder2;
        }

        PathResult Generate(double T, double dt, double S, double K, int lookback = 0) //create the path starting from S0, save the path as the array member in a new instance of the class PathResult, if lookback is 1 create a path for look back option pricing
        {
            // exception handling
            if (T <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Values T, K and S need to be greater than 0.");

            // exception handling
            if (lookback != 0 && lookback != 1)
                throw new System.ArgumentException("Value lookback takes values 0 and 1 only.");

            int m = Convert.ToInt32(T * (1 / dt)); //number of time steps

            //initialise arrays for storing output
            double[] Xn = new double[m + 1];
            double[] Sn = new double[m + 1];
            Xn[0] = Math.Log(S);
            Sn[0] = S;

            // Smin is changed only when calculating look back option price 
            double Smin = K;


            // create path
            for (int i = 0; i < m; i++)
            {
                Xn[i + 1] = Xn[i] + (r - 0.5 * SigmaDupSquared(T - i * dt, Math.Exp(Xn[i]), Smin)) * dt + Math.Sqrt(SigmaDupSquared(T - i * dt, Math.Exp(Xn[i]), Smin) * dt) * Normal.Sample(0, 1);
                Sn[i + 1] = Math.Exp(Xn[i + 1]);

                // update Smin if pricing lookback options
                if (Sn[i + 1] < Smin & lookback == 1) Smin = Sn[i + 1];
            }

            return new PathResult { Array = Sn };
        }

        public double MonteCarlo_EuropeanCallOptionPrice(double T, int N, double dt, double S, double K) //Monte Carlo method for pricing European Call Options
        {
            // exception handling
            if (T <= 0 || N <= 0 || dt <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("All values must be greater than 0.");

            double total = 0;
            //perform N iterations
            for (int i = 0; i < N; i++)
            {
                PathResult path = Generate(T, dt, S, K);

                total += Math.Exp(-r * T) * Math.Max(path.Array[path.Array.Length - 1] - K, 0);
            }

            double cN = total * (1 / Convert.ToDouble(N));

            return cN;

        }
        public double MonteCarlo_EuropeanPutOptionPrice(double T, int N, double dt, double S, double K) // Monte Carlo method for calculating European put option price from the call option price
        {
            // exception handling
            if (T <= 0 || N <= 0 || dt <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("All values must be greater than 0.");

            double total = 0;
            //perform N iterations
            for (int i = 0; i < N; i++)
            {
                PathResult path = Generate(T, dt, S, K);

                total += Math.Exp(-r * T) * Math.Max(K - path.Array[path.Array.Length - 1], 0);
            }

            double pN = total * (1 / Convert.ToDouble(N));

            return pN;
        }

        public double AsianArithmeticCallOptionPrice(double T, int N, double dt, double[] Tvals, double S, double K) //Monte Carlo method for calculating Asian arithmetic call option price
        {
            // exception handling
            if (T <= 0 || N <= 0 || dt <= 0 || K <= 0 || S <= 0 || Tvals.Any(x => x <= 0))
                throw new System.ArgumentException("All values must be greater than 0.");

            // exception handling
            if (Tvals.Length == 0)
                throw new System.ArgumentException("Tvals must contain at least one time.");

            double total = 0;
            // perform N iterations
            for (int i = 0; i < N; i++)
            {
                PathResult path = Generate(T, dt, S, K);
                double sum = 0;
                double avg;
                //calculate the average of all Si obtained at each Ti value
                for (int j = 0; j < Tvals.Length; j++)
                {
                    // get the index closest to the Ti value, path may not pass through each Ti value exactly
                    int k = Convert.ToInt32(Tvals[j] / dt);
                    sum += path.Array[k];
                }
                avg = sum / Convert.ToDouble(Tvals.Length);

                //total the payoffs
                total += Math.Exp(-r * T) * Math.Max(avg - K, 0);
            }
            //find the average
            double cN = total * (1 / Convert.ToDouble(N));

            return cN;
        }

        public double AsianArithmeticPutOptionPrice(double T, int N, double timestep, double[] T_1toM, double S, double K) //Monte Carlo method for calculating Asian arithmetic put option price
        {
            // exception handling
            if (T <= 0 || N <= 0 || timestep <= 0 || K <= 0 || S <= 0 || T_1toM.Any(x => x <= 0))
                throw new System.ArgumentException("Need T > 0, N > 0, timestep > 0, K > 0 and S > 0.");

            // exception handling
            if (T_1toM.Length == 0)
                throw new System.ArgumentException("Tvals must contain at least one time.");

            double total = 0;
            // iterate N times
            for (int i = 0; i < N; i++)
            {
                PathResult path = Generate(T, timestep, S, K);

                double sum = 0;
                double avg;
                //calculate the average of all Si obtained at each Ti value
                for (int j = 0; j < T_1toM.Length; j++)
                {
                    // get the index closest to the Ti value, path may not pass through each Ti value exactly
                    int k = Convert.ToInt32(T_1toM[i] / timestep);
                    sum += path.Array[k];
                }

                avg = sum / Convert.ToDouble(T_1toM.Length);

                //total the payoffs
                total += Math.Exp(-r * T) * Math.Max(K - avg, 0);
            }

            //find the average
            double pN = total * (1 / Convert.ToDouble(N));

            return pN;
        }

        public double LookBackCallOptionPrice(double T, int N, double timestep, double S, double K) //monte carlo method for finding the look back call option price
        {
            // exception handling
            if (T <= 0 || N <= 0 || timestep <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Values must be greater than 0.");

            double total = 0;
            
            //perform N iterations
            for (int i = 0; i < N; i++)
            {
                PathResult path = Generate(T, timestep, S, K, 1);

                // find the minimun value in the path
                double min = path.Array[0];
                foreach (double value in path.Array)
                {
                    if (value < min) min = value;
                }

                //total the payoffs
                total += Math.Exp(-r * T) * (path.Array[path.Array.Length - 1] - min);
            }

            //find the average
            double cN = total * (1 / Convert.ToDouble(N));

            return cN;
        }

        //class for storing the result of running the path method
        public class PathResult
        {
            public double[] Array { get; set;}
        }

    }
        
    class Program
    {
        static void Main(string[] args)
        {
            BlackScholesFormula task1 = new BlackScholesFormula();
            BlackScholesMonteCarlo MC = new BlackScholesMonteCarlo();

            double c;
            double p;
            double S = 100;
            int N = 100;
            double dt = 0.0078125;

            double[] K_table1 = new double[3] { 102.5313, 105.1271, 103.8212 };
            double[] T_table1 = new double[3] { 1, 2, 1.5 };
            double[] K_table2 = new double[5] { 80, 80, 90, 120, 90 };
            double[] T_table2 = new double[5] { 1, 2, 1, 1, 2 };

            //Task 2.1
            //calculate prices in Table 1
            double price_table1_1_task1 = task1.CalculateCallOptionPrice_atmoney(S, K_table1[0], T_table1[0]);
            double price_table1_2_task1 = task1.CalculateCallOptionPrice_atmoney(S, K_table1[1], T_table1[1]);
            double price_table1_3_task1 = task1.CalculateCallOptionPrice_atmoney(S, K_table1[2], T_table1[2]);

            //calculate prices in Table 2
            double price_table2_1_task1 = task1.CalculateCallOptionPrice(S, K_table2[0], T_table2[0]);
            double price_table2_2_task1 = task1.CalculateCallOptionPrice(S, K_table2[1], T_table2[1]);
            double price_table2_3_task1 = task1.CalculateCallOptionPrice(S, K_table2[2], T_table2[2]);
            double price_table2_4_task1 = task1.CalculatePutOptionPrice(S, K_table2[3], T_table2[3]);
            double price_table2_5_task1 = task1.CalculatePutOptionPrice(S, K_table2[4], T_table2[4]);

            //output prices in task 2.1
            Console.WriteLine("European Call Option price (when K=102.5313 and T=1)   = {0:f2}", price_table1_1_task1);
            Console.WriteLine("European Call Option price (when K=105.1271 and T=2)   = {0:f2}", price_table1_2_task1);
            Console.WriteLine("European Call Option price (when K=103.8212 and T=1.5) = {0:f2}", price_table1_3_task1);

            Console.WriteLine("European Call Option price (when K=80 and T=1)  = {0:f2}", price_table2_1_task1);
            Console.WriteLine("European Call Option price (when K=80 and T=2)  = {0:f2}", price_table2_2_task1);
            Console.WriteLine("European Call Option price (when K=90 and T=1)  = {0:f2}", price_table2_3_task1);
            Console.WriteLine("European Put Option price  (when K=120 and T=1) = {0:f2}", price_table2_4_task1);
            Console.WriteLine("European Put Option price  (when K=90 and T=2)  = {0:f2}", price_table2_5_task1);


            //Task 2.2
            //calculate prices in Table 1
            double price_table1_1_task2 = MC.MonteCarlo_EuropeanCallOptionPrice(T_table1[0], N, dt, S, K_table1[0]);
            double price_table1_2_task2 = MC.MonteCarlo_EuropeanCallOptionPrice(T_table1[0], N, dt, S, K_table1[1]);
            double price_table1_3_task2 = MC.MonteCarlo_EuropeanCallOptionPrice(T_table1[0], N, dt, S, K_table1[2]);

            //calculate prices in Table 2
            double price_table2_1_task2 = MC.MonteCarlo_EuropeanCallOptionPrice(T_table1[0], N, dt, S, K_table1[0]);
            double price_table2_2_task2 = MC.MonteCarlo_EuropeanCallOptionPrice(T_table1[0], N, dt, S, K_table1[1]);
            double price_table2_3_task2 = MC.MonteCarlo_EuropeanCallOptionPrice(T_table1[0], N, dt, S, K_table1[2]);
            double price_table2_4_task2 = MC.MonteCarlo_EuropeanPutOptionPrice(T_table1[0], N, dt, S, K_table1[3]);
            double price_table2_5_task2 = MC.MonteCarlo_EuropeanPutOptionPrice(T_table1[0], N, dt, S, K_table1[0]);

            //output prices in task 2.2
            Console.WriteLine("European Call Option price (Monte Carlo) (when K=102.5313 and T=1)   = {0:f2}", price_table1_1_task2);
            Console.WriteLine("European Call Option price (Monte Carlo) (when K=105.1271 and T=2)   = {0:f2}", price_table1_2_task2);
            Console.WriteLine("European Call Option price (Monte Carlo) (when K=103.8212 and T=1.5) = {0:f2}", price_table1_3_task2);

            Console.WriteLine("European Call Option price (Monte Carlo) (when K=80 and T=1)  = {0:f2}", price_table2_1_task2);
            Console.WriteLine("European Call Option price (Monte Carlo) (when K=80 and T=2)  = {0:f2}", price_table2_2_task2);
            Console.WriteLine("European Call Option price (Monte Carlo) (when K=90 and T=1)  = {0:f2}", price_table2_3_task2);
            Console.WriteLine("European Put Option price  (Monte Carlo) (when K=120 and T=1) = {0:f2}", price_table2_4_task2);
            Console.WriteLine("European Put Option price  (Monte Carlo) (when K=90 and T=2)  = {0:f2}", price_table2_5_task2);


            //Task 2.4
            double[] Tvals1 = new double[5] { 0.1, 0.2, 0.3, 0.4, 0.5 };
            double[] Tvals2 = new double[4] { 0.25, 0.5, 0.75, 1 };
            double[] Tvals3 = new double[2] { 0.5, 1 };

            Console.WriteLine("Asian Arithmetic Call Options");
            Console.WriteLine("Task 4 Table 1");
            Console.WriteLine("{0,-10}||{1,10}||{2,10}", "Strike", "Exercise", "Price");

            try
            {
                c = MC.AsianArithmeticCallOptionPrice(0.5, N, dt, Tvals1, S, 100);

                Console.WriteLine("{0,-10}||{1,10}||{2,10:n}", 100, 0.5, c);

                c = MC.AsianArithmeticCallOptionPrice(1, N, dt, Tvals2, S, 100);

                Console.WriteLine("{0,-10}||{1,10}||{2,10:n}", 100, 1, c);

                c = MC.AsianArithmeticCallOptionPrice(1.5, N, dt, Tvals3, S, 100);

                Console.WriteLine("{0,-10}||{1,10}||{2,10:n}", 100, 1.5, c);
            }
            catch (Exception e)
            {
                System.Console.WriteLine("Error: " + e.Message);
            }


            // Task 2.5
            Console.WriteLine("Look Back Call Options Monte Carlo with Floating Strike");
            Console.WriteLine("Task 5 Table 1");
            Console.WriteLine("{0,-10}||{1,10}", "Exercise", "Price");

            double[] Tval3 = new double[3] { 0.5, 1, 1.5 };

            try
            {
                for (int i = 0; i < 3; i++)
                {
                    c = MC.LookBackCallOptionPrice(Tval3[i], N, dt, S, 100);
                    Console.WriteLine("{0,-10}||{1,10:n}", Tval3[i], c);
                }
            }
            catch (Exception e)
            {
                System.Console.WriteLine("Error: " + e.Message);
            }

            Console.ReadKey();
        }
    }
}