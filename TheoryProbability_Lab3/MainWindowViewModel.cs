using LiveCharts;
using LiveCharts.Defaults;
using LiveCharts.Wpf;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using mathnet = MathNet.Numerics;

namespace TheoryProbability_Lab3
{
    public class MainWindowViewModel : ViewModel
    {
        public Func<double, string> YFormatter { get; set; }
        private static List<double>? _frequencies;
        private static List<double>? _randoms;
        private static string? _samplingCount;
        private static string? _paramB;
        private static string? _paramA;
        private static string? _beansCount;
        private static string? _h;
        private static SeriesCollection? _mainSeriesCollectionList;
        private static SeriesCollection? _optionalSeriesCollectionList;
        private static int? _distributionSelectedIndex;
        private static int? _kernelSelectedIndex;
        private static bool _isKdEPlotBuild;
        private static bool _isPointDensityBuild;
        private static bool _isSecondParamAvailable;
        private static double _maxValue;
        private static Visibility _isOptionalGraphVisible;
        private static List<string> _samplingInfoList;

        private const double PADDING = 0.02;

        public ICommand CreateDistributionGraphCommand { get; set; } // Команда обработки действия по нажатию на кнопку построения графика по новой выборке
        public ICommand CreateDistributionGraphOnCurrentCommand { get; set; } // Команда обработки действия по нажатию на кнопку построения графика по текущей выборке
        public string? SamplingCount // Объем выборки
        {
            get => _samplingCount;
            set => Set(ref _samplingCount, value);
        }
        public string? ParamA // Для нормального распределения - мат.ожидание, для равномерного - начало диапазона, для экспоненциального - параметр не учитывается
        {
            get => _paramA;
            set => Set(ref _paramA, value);
        }
        public string? ParamB // Для нормального распределения - дисперсия, для равномерного - конец диапазона, для экспоненциального - лямбда
        {
            get => _paramB;
            set => Set(ref _paramB, value);
        }
        public string? BeansCount // Количество интервалов
        {
            get => _beansCount;
            set => Set(ref _beansCount, value);
        }
        public string? H // Параметр размытости
        {
            get => _h;
            set => Set(ref _h, value);
        }
        public double MaxValue // Максимальное значение диапазона по оси ординат
        {
            get => _maxValue;
            set => Set(ref _maxValue, value);
        }
        public bool IsKdEPlotBuild // Переменная, отвечающая за построение графика оценки плотности распределения
        {
            get => _isKdEPlotBuild;
            set => Set(ref _isKdEPlotBuild, value);
        }
        public bool IsPointDensityBuild // Переменная, отвечающая за построение графика плотности точек
        {
            get => _isPointDensityBuild;
            set
            {
                IsOptionalGraphVisible = value ? Visibility.Visible : Visibility.Collapsed;
                Set(ref _isPointDensityBuild, value);
            }
        }
        public bool IsSecondParamAvailable // Переменная, отвечающая за доступность параметра A при экспоненциальном распределении
        {
            get => _isSecondParamAvailable;
            set => Set(ref _isSecondParamAvailable, value);
        }
        public int? DistributionSelectedIndex // Выбранный индекс в списке распределений
        {
            get => _distributionSelectedIndex;
            set
            {
                IsSecondParamAvailable = value != 2;
                Set(ref _distributionSelectedIndex, value);
            }
        }
        public int? KernelSelectedIndex // Выбранный индекс в списке ядерных функций
        {
            get => _kernelSelectedIndex;
            set => Set(ref _kernelSelectedIndex, value);
        }
        public SeriesCollection? MainSeriesCollectionList // Коллекция серий в LiveChart
        {
            get => _mainSeriesCollectionList;
            set => Set(ref _mainSeriesCollectionList, value);
        }
        public SeriesCollection? OptionalSeriesCollectionList // Коллекция серий в LiveChart
        {
            get => _optionalSeriesCollectionList;
            set => Set(ref _optionalSeriesCollectionList, value);
        }
        public Visibility IsOptionalGraphVisible
        {
            get => _isOptionalGraphVisible;
            set => Set(ref _isOptionalGraphVisible, value);
        }
        public List<string> SamplingInfoList
        {
            get => _samplingInfoList;
            set => Set(ref _samplingInfoList, value);
        } 

        private void CreateDistributionGrapOnCurrent_Executed(object obj)
        {
            if (SamplingCount == null || ParamA == null || ParamB == null || BeansCount == null || H == null) return;

            MainSeriesCollectionList?[2].Values.Clear();

            switch (DistributionSelectedIndex)
            {
                case 0: // Нормальное распределение
                    CreateNormalDistributionGraph(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount),
                        double.Parse(H.Replace(",","."), NumberStyles.Any, CultureInfo.InvariantCulture), true);
                    break;

                case 1: // Равномерное распределение
                    CreateUniformDistributionGraph(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount),
                        double.Parse(H.Replace(",", "."), NumberStyles.Any, CultureInfo.InvariantCulture), true);
                    break;

                case 2: // Экспоненциальное распределение
                    if (double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture) <= 0) MessageBox.Show("Введите положительное число в качестве параметра!");
                    else
                    {
                        CreateExponentialDistributionGraph(int.Parse(SamplingCount),
                            double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount),
                            double.Parse(H.Replace(",", "."), NumberStyles.Any, CultureInfo.InvariantCulture), true);
                    }
                    break;
            }
        }

        private void CreateDistributionGraph_Executed(object obj)
        {
            if (SamplingCount == null || ParamA == null || ParamB == null || BeansCount == null || H == null) return;

            MainSeriesCollectionList?[0].Values.Clear();
            MainSeriesCollectionList?[1].Values.Clear();
            MainSeriesCollectionList?[2].Values.Clear();

            OptionalSeriesCollectionList?[0].Values.Clear();
            switch (DistributionSelectedIndex)
            {
                case 0: // Нормальное распределение
                    CreateNormalDistributionGraph(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount), 
                        double.Parse(H.Replace(",", "."), NumberStyles.Any, CultureInfo.InvariantCulture), false);
                    break;

                case 1: // Равномерное распределение
                    CreateUniformDistributionGraph(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount),
                        double.Parse(H.Replace(",", "."), NumberStyles.Any, CultureInfo.InvariantCulture), false);
                    break;

                case 2: // Экспоненциальное распределение
                    if (double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture) <= 0) MessageBox.Show("Введите положительное число в качестве параметра!");
                    else
                    {
                        CreateExponentialDistributionGraph(int.Parse(SamplingCount),
                            double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount),
                            double.Parse(H.Replace(",", "."), NumberStyles.Any, CultureInfo.InvariantCulture), false);
                    }
                    break;
            }
        }
        public MainWindowViewModel()
        {
            YFormatter = value => value.ToString("N");
            MainSeriesCollectionList = new SeriesCollection();
            OptionalSeriesCollectionList = new SeriesCollection();
            SamplingCount = "100";
            ParamA = "0";
            ParamB = "1";
            MaxValue = 1.0;
            BeansCount = "10";
            H = "1";
            DistributionSelectedIndex = 0;
            KernelSelectedIndex = 0;
            IsSecondParamAvailable = true;
            IsOptionalGraphVisible = Visibility.Collapsed;

            CreateDistributionGraphCommand = new Command(CreateDistributionGraph_Executed, _ => true);
            CreateDistributionGraphOnCurrentCommand = new Command(CreateDistributionGrapOnCurrent_Executed, _ => true);

            MainSeriesCollectionList.Add(new ColumnSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                MaxColumnWidth = 45,
                Stroke = new SolidColorBrush(Color.FromRgb(10, 15, 100)),
                StrokeThickness = 1.5
            });

            MainSeriesCollectionList.Add(new LineSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                StrokeThickness = 3,
            });

            MainSeriesCollectionList.Add(new LineSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                StrokeThickness = 3,
                Stroke = new SolidColorBrush(Colors.DarkSlateBlue),
                Fill = new SolidColorBrush(Color.FromArgb(30, 10, 54, 166))
            });

            OptionalSeriesCollectionList.Add(new LineSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                StrokeThickness = 3,
                Stroke = new SolidColorBrush(Colors.ForestGreen),
                PointGeometry = DefaultGeometries.Diamond,
                Fill = new SolidColorBrush(Colors.Transparent)

            });
        }

        private void GetSamplingInfo(double hi, double moda, double maxPvalue, double degOfFreedom)
        {
            if (_randoms == null) return;

            SamplingInfoList = new List<string>();
            var randoms = _randoms;
            randoms.Sort();
            var median = MathNet.Numerics.Statistics.Statistics.Median(randoms);

            var sampleAverage = randoms.Sum()/randoms.Count;
            var stDev = MathNet.Numerics.Statistics.Statistics.StandardDeviation(randoms);

            var range = randoms.Max() - randoms.Min();

            var lowerConfBound = sampleAverage - 1.96 * (stDev / Math.Sqrt(randoms.Count));
            var upperConfBound = sampleAverage + 1.96 * stDev / Math.Sqrt(randoms.Count);


            SamplingInfoList = new List<string>
            {
                $"Выборочное среднее: {Math.Round(sampleAverage, 4)}",
                $"Медиана: {Math.Round(median, 4)}",
                $"Дисперсия: {Math.Round(stDev, 4)}",
                $"Размах: {Math.Round(range, 4)}",
                $"Мода: {Math.Round(moda, 4)}",
                $"Довер. интервал: {Math.Round(lowerConfBound, 4)} ± {Math.Round(upperConfBound, 4)}",
                $"Хи-квадрат: {Math.Round(hi, 4)}",
                $"Критическая точка хи-распределения при уровне значимости 0,95 и кол-ве степеней свободы = {degOfFreedom}, равна: {Math.Round(maxPvalue,4)}"
            };
            
        }

        #region Нормальное распределение
        private void CreateNormalDistributionGraph(int samplingCount, double mean, double stDev, int beansCount, double h, bool inCurrent) // Построение графика нормального распределения
        {
            var min3Sigm = mean - 3 * stDev;
            var max3Sigm = mean + 3 * stDev;
            var incr = Math.Abs(max3Sigm - min3Sigm) / beansCount;
            if (inCurrent)
            {
                if(IsKdEPlotBuild)
                    KernelDensityGraphBuild(min3Sigm, max3Sigm, incr, h, samplingCount);
                return;
            }

            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++) 
            {
                _randoms[i] = Calculus.GenerateNormalDistributionRandomValue(mean, stDev);
            }

            var yMax = 0.0;

            var vals = new List<double>();
            var hiSquare = 0.0;
            var moda = 0.0;
            var maxFreq = 0.0;
            for (var i = min3Sigm; i <= max3Sigm; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(c => c >= i && c <= i + incr) / (double)(samplingCount);
                _frequencies[valsCount] = count;

                MainSeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));

                var normalDensityRandomValue = incr * Calculus.NormalDistributionDensityFunction(stDev, mean, i);
                
                MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(i, normalDensityRandomValue));

                if (valsCount == beansCount - 1)
                    MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(i + incr, incr * Calculus.NormalDistributionDensityFunction(stDev, mean, i + incr)));

                yMax = normalDensityRandomValue > yMax ? normalDensityRandomValue : yMax;
                vals.Add(i);

                hiSquare += Math.Pow(incr*Calculus.NormalDistributionDensityFunction(stDev, mean, i) - _frequencies[valsCount], 2) /
                   incr* Calculus.NormalDistributionDensityFunction(stDev, mean, i);

                if (!(_frequencies[valsCount] > maxFreq)) continue;
                maxFreq = _frequencies[valsCount];
                moda = i;
            }
            var maxPvalue = mathnet.Distributions.ChiSquared.InvCDF(beansCount-3, 0.05);
            GetSamplingInfo(hiSquare, moda, maxPvalue, beansCount-3);


            MaxValue = yMax > _frequencies.Max() ? yMax + PADDING : _frequencies.Max() + PADDING;
            if (IsPointDensityBuild)
            {
                _randoms.ForEach(rand => OptionalSeriesCollectionList?[0].Values.Add(new ObservablePoint(rand, 0)));
            }


            if (!IsKdEPlotBuild) return;
            KernelDensityGraphBuild(min3Sigm, max3Sigm, incr, h, samplingCount);


        }

        #endregion


        #region Равномерное распределение
        
        private void CreateUniformDistributionGraph(int samplingCount, double a, double b, int beansCount, double h, bool inCurrent) // Построение графика
        {
            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++)
            {
                _randoms[i] = Calculus.GenerateUniformDistributionRandomValue(a, b);
            }

            _randoms.Sort();
            var incr = Math.Abs(b - a) / beansCount;
            var yValue = 1 / (b - a);
            var vals = new List<double>();
            var yMax = 0.0;
            var maxFreq = 0.0;
            var moda = 0.0;
            var hiSquare = 0.0;
            for (var i = a; i <= b; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(rand => rand >= i && rand <= i + incr) / (double)samplingCount / incr;
                _frequencies[valsCount] = count;

                MainSeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));
                vals.Add(i);
                MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(i, yValue));
                yMax = yValue > yMax ? yValue : yMax;
                hiSquare += Math.Pow(yValue -_frequencies[valsCount], 2) / yValue;
                if (!(_frequencies[valsCount] > maxFreq)) continue;
                maxFreq = _frequencies[valsCount];
                moda = i;
            }
            
            MaxValue = yMax > _frequencies.Max() ? yMax + PADDING : _frequencies.Max() + PADDING;

            if (IsPointDensityBuild)
            {
                _randoms.ForEach(rand => OptionalSeriesCollectionList?[0].Values.Add(new ObservablePoint(rand, 0)));
            }
            var maxPvalue = mathnet.Distributions.ChiSquared.InvCDF(beansCount - 3, 0.05);
            GetSamplingInfo(hiSquare, moda, maxPvalue, beansCount-3);

            if (!IsKdEPlotBuild) return;

            KernelDensityGraphBuild(a, b, incr, h, samplingCount);

        }


        #endregion


        #region Экспоненциальное распределение
        
        private void CreateExponentialDistributionGraph(int samplingCount, double lambda, int beansCount, double h, bool inCurrent) // Построение графика
        {
            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++)
            {
                _randoms[i] = Calculus.GenerateExponentialDistributionRandomValue(lambda);
            }

            _randoms.Sort();
            var minVal = _randoms[0];
            var maxVal = _randoms[^1];
            var incr = Math.Abs(maxVal - minVal) / beansCount;
            var maxFreq = 0.0;
            var moda = 0.0;
            var hiSquare = 0.0;
            var vals = new List<double>();
            var yMax = 0.0;
            for (var i = minVal; i <= maxVal; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(c => c >= i && c <= i + incr) / (double)samplingCount;
                _frequencies[valsCount] = count;
                MainSeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));
                vals.Add(i);
                var yValue = Calculus.ExponentialDistributionDensityFunction(lambda, i);

                MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(i, yValue));
                yMax = yValue > yMax ? yValue : yMax;
                hiSquare += Math.Pow(Calculus.ExponentialDistributionDensityFunction(lambda, i) - 
                    _frequencies[valsCount], 2) / Calculus.ExponentialDistributionDensityFunction(lambda, i);
                if (!(_frequencies[valsCount] > maxFreq)) continue;
                maxFreq = _frequencies[valsCount];
                moda = i;
            }
            
            MaxValue = yMax > _frequencies.Max() ? yMax + PADDING : _frequencies.Max() + PADDING;
            var maxPvalue = mathnet.Distributions.ChiSquared.InvCDF(beansCount - 3, 0.05);
            GetSamplingInfo(hiSquare, moda, maxPvalue, beansCount - 3);
            if (IsPointDensityBuild)
            {
                _randoms.ForEach(rand => OptionalSeriesCollectionList?[0].Values.Add(new ObservablePoint(rand, 0)));
            }
            if (!IsKdEPlotBuild) return;

            KernelDensityGraphBuild(minVal, maxVal, incr, h, samplingCount);
       
        }
        
        #endregion
        private void KernelDensityGraphBuild(double min, double max, double incr, double h, double samplingCount)
        {
            if (_randoms == null) return;
            var upperBound = DistributionSelectedIndex == 0 ? max + incr : max;
            var multiplyKoeff = DistributionSelectedIndex == 0 ? incr : 1.0;
            switch (KernelSelectedIndex)
            {
                case 0:
                    for (var i = min; i <= upperBound; i += incr)
                    {
                        var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => Calculus.GaussianKernel((i - rand) / h));
                        MainSeriesCollectionList?[2].Values.Add(new ObservablePoint(i, multiplyKoeff * kde));
                    }
                    break;

                case 1:
                    for (var i = min; i <= upperBound; i += incr)
                    {
                        var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => Calculus.EpanechnikovKernel((i - rand) / h));
                        MainSeriesCollectionList?[2].Values.Add(new ObservablePoint(i, multiplyKoeff * kde));
                    }
                    break;

                case 2:
                    for (var i = min; i <= upperBound; i += incr)
                    {
                        var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => Calculus.LogisticKernel((i - rand) / h));
                        MainSeriesCollectionList?[2].Values.Add(new ObservablePoint(i, multiplyKoeff * kde));
                    }
                    break;
            }
        }
        
    }
}
