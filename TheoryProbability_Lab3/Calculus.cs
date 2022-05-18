using System;
using static System.Math;
namespace TheoryProbability_Lab3
{
    internal static class Calculus
    {
        private static readonly Random Rand = new();

        #region Генерация чисел
        public static double GenerateNormalDistributionRandomValue(double mean, double stDev) // Генерация случайного числа по нормальному закону распределения
        {
            var u1 = Rand.NextDouble();
            var u2 = Rand.NextDouble();
            var randStdNormal = Sqrt(-2.0 * Log(u1)) * Sin(2.0 * PI * u2);
            var randNormal = mean + stDev * randStdNormal;
            return randNormal;
        }

        public static double GenerateUniformDistributionRandomValue(double start, double stop) // Генерация случайного числа по равномерному закону распределения
        {
            return start + Rand.NextDouble() * (stop - start);
        }
        public static double GenerateExponentialDistributionRandomValue(double lambda) // Генерация случайного числа по экспоненциальному закону распределения
        {
            var uniformValue = Rand.NextDouble();
            return (-1 / lambda) * Log(1 - uniformValue);
        }

        #endregion


        #region Плотность распределений
        public static double NormalDistributionDensityFunction(double sigm, double mean, double x) // Плотность нормального распределния
        {
            return (1 / (sigm * Sqrt(2 * PI))) * Pow(E, (-1 * ((Pow(x - mean, 2)) / (2 * sigm * sigm))));
        }
        public static double ExponentialDistributionDensityFunction(double lambda, double x) // Плотность экспоненциального распределения
        {
            return lambda * Exp(-lambda * x);
        }

        #endregion


        #region Ядерные функции
        public static double GaussianKernel(double x) // Гауссово ядро
        {
            return Exp((-1 * x * x) / 2) / Sqrt(2 * PI);
        }
        public static double EpanechnikovKernel(double x) // Епанечниково ядро
        {
            return Abs(x) <= 1 ? 0.75 * (1 - x * x) : 0.0;
        }
        public static double LogisticKernel(double x) // Логистическое ядро
        {
            return 1 / (Exp(x) + 2 + Exp(-x));
        }


        #endregion

    }
}
