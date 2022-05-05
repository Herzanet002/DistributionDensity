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

namespace TheoryProbability_Lab3
{
    public class MainWindowViewModel : ViewModel
    {
        public Func<double, string> YFormatter { get; set; }
        private readonly Random _rand = new();
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

        private const double PADDING = 0.02;

        public ICommand CreateNormalDistributionPlotCommand { get; set; } // Команда обработки действия по нажатию на кнопку построения графика по новой выборке
        //public ICommand CreateNormalDistributionPlotCommand { get; set; } // Команда обработки действия по нажатию на кнопку построения графика по текущей выборке
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
        private void CreateDistributionPlot_Executed(object obj)
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
                        double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount), double.Parse(H, NumberStyles.Any, CultureInfo.InvariantCulture));
                    break;

                case 1: // Равномерное распределение
                    CreateUniformDistributionGraph(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount), double.Parse(H, NumberStyles.Any, CultureInfo.InvariantCulture));
                    break;

                case 2: // Экспоненциальное распределение
                    if (double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture) <= 0) MessageBox.Show("Введите положительное число в качестве параметра!");
                    else
                    {
                        CreateExponentialDistributionGraph(int.Parse(SamplingCount),
                            double.Parse(ParamB, NumberStyles.Any, CultureInfo.InvariantCulture), int.Parse(BeansCount), double.Parse(H, NumberStyles.Any, CultureInfo.InvariantCulture));
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

            CreateNormalDistributionPlotCommand = new Command(CreateDistributionPlot_Executed, _ => true);

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


        #region Нормальное распределение
        private double GenerateNormalDistributionRandomValue(double mean, double stDev) // Генерация случайного числа по нормальному закону распределения
        {
            var u1 = _rand.NextDouble();
            var u2 = _rand.NextDouble();
            var randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
            var randNormal = mean + stDev * randStdNormal;
            return randNormal;
        }
        private void CreateNormalDistributionGraph(int samplingCount, double mean, double stDev, int beansCount, double h) // Построение графика
        {
            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++) //выборка по нормальному
            {
                _randoms[i] = GenerateNormalDistributionRandomValue(mean, stDev);
            }
            var yMax = 0.0;
            // _randoms.Sort();
            var min3Sigm = mean - 3 * stDev;
            var max3Sigm = mean + 3 * stDev;
            var incr = Math.Abs(max3Sigm - min3Sigm) / beansCount;

            var vals = new List<double>();

            for (var i = min3Sigm; i <= max3Sigm; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(c => c >= i && c <= i + incr) / (double)(samplingCount);
                _frequencies[valsCount] = count;

                MainSeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));

                var normalDensityRandomValue = incr * NormalDistributionDensityFunction(stDev, mean, i);

                MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(i, normalDensityRandomValue));

                if (valsCount == beansCount - 1)
                    MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(i + incr, incr * NormalDistributionDensityFunction(stDev, mean, i + incr)));

                yMax = normalDensityRandomValue > yMax ? normalDensityRandomValue : yMax;
                vals.Add(i);
            }



            MaxValue = yMax > _frequencies.Max() ? yMax + PADDING : _frequencies.Max() + PADDING;
            if (IsPointDensityBuild)
            {
                _randoms.ForEach(rand => OptionalSeriesCollectionList?[0].Values.Add(new ObservablePoint(rand, 0)));

            }


            if (!IsKdEPlotBuild) return;
            KernelDensityGraphBuild(min3Sigm, max3Sigm, incr, h, samplingCount);


        }


        private static double NormalDistributionDensityFunction(double sigm, double mean, double x) // Плотность нормального распределния
        {
            return (1 / (sigm * Math.Sqrt(2 * Math.PI))) * Math.Pow(Math.E, (-1 * ((Math.Pow(x - mean, 2)) / (2 * sigm * sigm))));
        }

        #endregion


        #region Равномерное распределение
        private double GenerateUniformDistributionRandomValue(double start, double stop) // Генерация случайного числа по равномерному закону распределения
        {
            return start + _rand.NextDouble() * (stop - start);
        }
        private void CreateUniformDistributionGraph(int samplingCount, double a, double b, int beansCount, double h) // Построение графика
        {
            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++)
            {
                _randoms[i] = GenerateUniformDistributionRandomValue(a, b);
            }

            _randoms.Sort();
            var minVal = a;
            var maxVal = b;
            var incr = Math.Abs(maxVal - minVal) / beansCount;

            var vals = new List<double>();
            for (var i = minVal; i <= maxVal; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(rand => rand >= i && rand <= i + incr) / (double)samplingCount / incr;
                _frequencies[valsCount] = count;

                MainSeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));
                vals.Add(i);

            }
            var yMax = 0.0;
            vals.ForEach(v =>
            {
                var yValue = 1 / (maxVal - minVal);
                MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(v, yValue));
                yMax = yValue > yMax ? yValue : yMax;
            });
            MaxValue = yMax > _frequencies.Max() ? yMax + PADDING : _frequencies.Max() + PADDING;

            if (IsPointDensityBuild)
            {
                _randoms.ForEach(rand => OptionalSeriesCollectionList?[0].Values.Add(new ObservablePoint(rand, 0)));
            }

            if (!IsKdEPlotBuild) return;

            KernelDensityGraphBuild(minVal, maxVal, incr, h, samplingCount);

        }


        #endregion


        #region Экспоненциальное распределение
        private double GenerateExponentialDistributionRandomValue(double lambda) // Генерация случайного числа по экспоненциальному закону распределения
        {
            var uniformValue = _rand.NextDouble();
            return (-1 / lambda) * Math.Log(1 - uniformValue);
        }
        private void CreateExponentialDistributionGraph(int samplingCount, double lambda, int beansCount, double h) // Построение графика
        {
            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++)
            {
                _randoms[i] = GenerateExponentialDistributionRandomValue(lambda);
            }

            _randoms.Sort();
            var minVal = _randoms[0];
            var maxVal = _randoms[^1];
            var incr = Math.Abs(maxVal - minVal) / beansCount;


            var vals = new List<double>();
            for (var i = minVal; i <= maxVal; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(c => c >= i && c <= i + incr) / (double)samplingCount;
                _frequencies[valsCount] = count;
                MainSeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));
                vals.Add(i);

            }
            var yMax = 0.0;

            vals.ForEach(v =>
            {
                var yValue = ExponentialDistributionDensityFunction(lambda, v);

                MainSeriesCollectionList?[1].Values.Add(new ObservablePoint(v, yValue));
                yMax = yValue > yMax ? yValue : yMax;
            });
            MaxValue = yMax > _frequencies.Max() ? yMax + PADDING : _frequencies.Max() + PADDING;
            if (IsPointDensityBuild)
            {
                _randoms.ForEach(rand => OptionalSeriesCollectionList?[0].Values.Add(new ObservablePoint(rand, 0)));
            }
            if (!IsKdEPlotBuild) return;

            KernelDensityGraphBuild(minVal, maxVal, incr, h, samplingCount);
       
        }
        private static double ExponentialDistributionDensityFunction(double lambda, double x) // Плотность жкспоненциального распределения
        {
            return lambda * Math.Exp(-lambda * x);
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
                        var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => GaussianKernel((i - rand) / h));
                        MainSeriesCollectionList?[2].Values.Add(new ObservablePoint(i, multiplyKoeff * kde));
                    }
                    break;

                case 1:
                    for (var i = min; i <= upperBound; i += incr)
                    {
                        var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => EpanechnikovKernel((i - rand) / h));
                        MainSeriesCollectionList?[2].Values.Add(new ObservablePoint(i, multiplyKoeff * kde));
                    }
                    break;

                case 2:
                    for (var i = min; i <= upperBound; i += incr)
                    {
                        var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => LogisticKernel((i - rand) / h));
                        MainSeriesCollectionList?[2].Values.Add(new ObservablePoint(i, multiplyKoeff * kde));
                    }
                    break;
            }
        }
        private static double GaussianKernel(double x) // Гауссово ядро
        {
            return Math.Exp((-1 * x * x) / 2) / Math.Sqrt(2 * Math.PI);
        }
        private static double EpanechnikovKernel(double x) // Епанечниково ядро
        {
            return Math.Abs(x) <= 1 ? 0.75 * (1 - x * x) : 0.0;
        }
        private static double LogisticKernel(double x) // Логистическое ядро
        {
            return 1 / (Math.Exp(x) + 2 + Math.Exp(-x));
        }
    }
}
