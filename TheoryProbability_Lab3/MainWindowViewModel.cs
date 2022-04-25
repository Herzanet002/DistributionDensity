using LiveCharts;
using LiveCharts.Defaults;
using LiveCharts.Wpf;
using System;
using System.Collections.Generic;
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
        private static SeriesCollection? _seriesCollectionList;
        private static int? _comboboxSelectedIndex;
        private static bool _isKdEPlotBuild;
        private static bool _isPointDensityBuild;
        private static bool _isSecondParamAvailable;
        private static double _maxValue;

        public ICommand CreateNormalDistributionPlotCommand { get; set; } // Команда обработки действия по нажатию на кнопку построения графика
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
            set => Set(ref _isPointDensityBuild, value);
        }

        public bool IsSecondParamAvailable // Переменная, отвечающая за доступность параметра A при экспоненциальном распределении
        {
            get => _isSecondParamAvailable;
            set => Set(ref _isSecondParamAvailable, value);
        }
        public int? ComboboxSelectedIndex // Выбранный индекс в списке распределений
        {
            get => _comboboxSelectedIndex;
            set
            {
                IsSecondParamAvailable = value != 2;
                Set(ref _comboboxSelectedIndex, value);
            }
        }

        public SeriesCollection? SeriesCollectionList // Коллекция серий в LiveChart
        {
            get => _seriesCollectionList;
            set => Set(ref _seriesCollectionList, value);
        }
        private void CreateDistributionPlot_Executed(object obj)
        {
            if (SamplingCount == null || ParamA == null || ParamB == null || BeansCount == null || H == null) return;

            SeriesCollectionList?[0].Values.Clear();
            SeriesCollectionList?[1].Values.Clear();
            SeriesCollectionList?[2].Values.Clear();
            SeriesCollectionList?[3].Values.Clear();
            switch (ComboboxSelectedIndex)
            {
                case 0: // Нормальное распределение
                    CreateNormalDistributionPlot(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB), int.Parse(BeansCount), double.Parse(H));
                    break;

                case 1: // Равномерное распределение
                    CreateUniformDistributionPlot(int.Parse(SamplingCount), double.Parse(ParamA),
                        double.Parse(ParamB), int.Parse(BeansCount), double.Parse(H));
                    break;

                case 2: // Экспоненциальное распределение
                    if (double.Parse(ParamB) <= 0) MessageBox.Show("Введите положительное число в качестве параметра!");
                    else
                    {
                        CreateExponentialDistributionPlot(int.Parse(SamplingCount), 
                            double.Parse(ParamB), int.Parse(BeansCount), double.Parse(H));
                    }
                    break;
            }
        }
        public MainWindowViewModel()
        {
            YFormatter = value => value.ToString("N");
            SeriesCollectionList = new SeriesCollection();
            SamplingCount = "100";
            ParamA = "0";
            ParamB = "1";
            MaxValue = 1.0;
            BeansCount = "10";
            H = "1";
            ComboboxSelectedIndex = 0;
            IsSecondParamAvailable = true;
            CreateNormalDistributionPlotCommand = new Command(CreateDistributionPlot_Executed, _ => true);
            SeriesCollectionList.Add(new ColumnSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                MaxColumnWidth = 45,
                Stroke = new SolidColorBrush(Color.FromRgb(10, 15, 100)),
                StrokeThickness = 1
            });
            SeriesCollectionList.Add(new LineSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                StrokeThickness = 3,

            });
            SeriesCollectionList.Add(new LineSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                StrokeThickness = 3
            });
            SeriesCollectionList.Add(new LineSeries
            {
                Values = new ChartValues<ObservablePoint>(),
                StrokeThickness = 3,
                Stroke = new SolidColorBrush(Colors.ForestGreen),
                PointGeometry = DefaultGeometries.Diamond,

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
        private void CreateNormalDistributionPlot(int samplingCount, double mean, double stDev, int beansCount, double h) // Построение графика
        {
            _frequencies = new List<double>(new double[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++) //выборка по нормальному
            {
                _randoms[i] = GenerateNormalDistributionRandomValue(mean, stDev);
            }
            var yMax = 0.0;
            _randoms.Sort();
            var minVal = mean - 3 * stDev;
            var maxVal = mean + 3 * stDev;
            var incr = Math.Abs(maxVal - minVal) / beansCount;

            var vals = new List<double>();
            for (var i = minVal; i <= maxVal; i += incr)
            {
                var valsCount = vals.Count;
                if (valsCount >= beansCount) break;
                var count = _randoms.Count(c => c >= i && c <= i + incr) / (double)(samplingCount);
                _frequencies[valsCount] = count;

                SeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));

                var y = incr * NormalDistributionDensityFunction(stDev, mean, i);

                SeriesCollectionList?[1].Values.Add(new ObservablePoint(i, y));

                if (valsCount == beansCount - 1)
                    SeriesCollectionList?[1].Values.Add(new ObservablePoint(i + incr, incr * NormalDistributionDensityFunction(stDev, mean, i + incr)));

                yMax = y > yMax ? y : yMax;
                vals.Add(i);
            }



            MaxValue = yMax > _frequencies.Max() ? yMax + 0.02 : _frequencies.Max() + 0.02;
            if (IsPointDensityBuild)
            {
                foreach (var rand in _randoms)
                {
                    SeriesCollectionList?[3].Values.Add(new ObservablePoint(rand, 0));
                }
            }


            if (!IsKdEPlotBuild) return;
            for (var i = minVal; i <= maxVal; i += incr)
            {
                var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => GaussianKernel((i - rand) / h));

                SeriesCollectionList?[2].Values.Add(new ObservablePoint(i, incr * kde));
            }

        }
        private static double NormalDistributionDensityFunction(double sigm, double mean, double x) // Плотность нормального распределния
        {
            return (1 / (sigm * Math.Sqrt(2 * Math.PI))) * Math.Pow(Math.E, (-1 * ((Math.Pow(x - mean, 2)) / (2 * sigm * sigm))));
        }

        #endregion


        #region Равномерное распределение
        private double GenerateUniformDistributionRandomValue(double a, double b) // Генерация случайного числа по равномерному закону распределения
        {
            return a + _rand.NextDouble() * (b - a);
        }
        private void CreateUniformDistributionPlot(int samplingCount, double a, double b, int beansCount, double h) // Построение графика
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
                var count = _randoms.Count(c => c >= i && c <= i + incr) / (double)samplingCount;
                _frequencies[valsCount] = count;

                SeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));
                vals.Add(i);

            }
            var yMax = 0.0;
            vals.ForEach(v =>
            {
                var y = 1 / (maxVal - minVal);
                SeriesCollectionList?[1].Values.Add(new ObservablePoint(v, y));
                yMax = y > yMax ? y : yMax;
            });
            MaxValue = yMax > _frequencies.Max() ? yMax + 0.05 : _frequencies.Max() + 0.05;

            if (IsPointDensityBuild)
            {
                _randoms.ForEach(r => SeriesCollectionList?[3].Values.Add(new ObservablePoint(r, 0)));
            }

            if (!IsKdEPlotBuild) return;

            for (var i = minVal; i <= maxVal; i += incr)
            {
                var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => GaussianKernel((i - rand) / h));
                SeriesCollectionList?[2].Values.Add(new ObservablePoint(i, kde));
            }

        }


        #endregion


        #region Экспоненциальное распределение
        private double GenerateExponentialDistributionRandomValue(double lambda) // Генерация случайного числа по экспоненциальному закону распределения
        {
            var u = _rand.NextDouble();
            return (-1 / lambda) * Math.Log(1 - u);
        }
        private void CreateExponentialDistributionPlot(int samplingCount, double lambda, int beansCount, double h) // Построение графика
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
                SeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[valsCount]));
                vals.Add(i);

            }
            var yMax = 0.0;

            vals.ForEach(v =>
            {
                //var yMean = (2 * v + incr) / 2;
                var y = ExponentialDistributionDensityFunction(lambda, v);

                SeriesCollectionList?[1].Values.Add(new ObservablePoint(v, y));
                yMax = y > yMax ? y : yMax;
            });
            MaxValue = yMax > _frequencies.Max() ? yMax + 0.05 : _frequencies.Max() + 0.05;
            if (IsPointDensityBuild)
            {
                foreach (var rand in _randoms)
                {
                    SeriesCollectionList?[3].Values.Add(new ObservablePoint(rand, 0));
                }
            }
            if (!IsKdEPlotBuild) return;

            for (var i = minVal; i <= maxVal; i += incr)
            {
                var kde = (1.0 / samplingCount / h) * _randoms.Sum(rand => GaussianKernel((i - rand) / h));

                SeriesCollectionList?[2].Values.Add(new ObservablePoint(i, kde * incr));
            }

        }
        private static double ExponentialDistributionDensityFunction(double lambda, double x) // Плотность жкспоненциального распределения
        {
            return lambda * Math.Exp(-lambda * x);
        }


        #endregion


        private static double GaussianKernel(double x) // Гауссово ядро
        {
            return Math.Exp((-1 * x * x) / 2) / Math.Sqrt(2 * Math.PI);
        }

    }
}
