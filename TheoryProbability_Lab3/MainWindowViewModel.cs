using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Input;
using LiveCharts;
using LiveCharts.Defaults;
using LiveCharts.Wpf;

namespace TheoryProbability_Lab3
{
    public class MainWindowViewModel : ViewModel
    {
        public Func<double, string> YFormatter { get; set; }
        private readonly Random _rand = new();
        private List<int> _frequencies;
        private List<double> _randoms;
        private static string? _samplingCount;
        private static string? _stDev;
        private static string? _mean;
        private static string? _beansCount;
        private static SeriesCollection? _seriesCollectionList;
        private static int? _comboboxSelectedIndex;
        public ICommand CreateNormalDistributionPlotCommand { get; set; }
        public string? SamplingCount
        {
            get => _samplingCount;
            set => Set(ref _samplingCount, value);
        }
        public string? StDev
        {
            get => _stDev;
            set => Set(ref _stDev, value);
        }
        public string? Mean
        {
            get => _mean;
            set => Set(ref _mean, value);
        }
        public string? BeansCount
        {
            get => _beansCount;
            set => Set(ref _beansCount, value);
        }

        public int? ComboboxSelectedIndex
        {
            get => _comboboxSelectedIndex;
            set => Set(ref _comboboxSelectedIndex, value);
        }
        public SeriesCollection? SeriesCollectionList
        {
            get => _seriesCollectionList;
            set => Set(ref _seriesCollectionList, value);
        }
        private void CreateDistributionPlot_Executed(object obj)
        {
            if (SamplingCount == null || Mean == null || StDev == null || BeansCount == null) return;

            SeriesCollectionList?[0].Values.Clear();
            SeriesCollectionList?[1].Values.Clear();
            switch (ComboboxSelectedIndex)
            {
                case 0: //нормальное распределение
                    CreateNormalDistributionPlot(int.Parse(SamplingCount), double.Parse(Mean),
                        double.Parse(StDev), int.Parse(BeansCount));
                    break;

                case 1: //равномерное распределение
                    CreateUniformDistributionPlot(int.Parse(SamplingCount), double.Parse(Mean),
                        double.Parse(StDev), int.Parse(BeansCount));
                    break;

                case 2: //экспоненциальное распределение
                    CreateExponentialDistributionPlot(int.Parse(SamplingCount), double.Parse(Mean),
                        double.Parse(StDev), int.Parse(BeansCount));
                    break;
            }
        }
        
#pragma warning disable CS8618
        public MainWindowViewModel()
        {
            YFormatter = value => value.ToString("N");
            SeriesCollectionList = new SeriesCollection();
            SamplingCount = "100";
            Mean = "0";
            StDev = "1";
            BeansCount = "10";
            ComboboxSelectedIndex = 0;
            CreateNormalDistributionPlotCommand = new Command(CreateDistributionPlot_Executed, _ => true);
            SeriesCollectionList.Add(new ColumnSeries
            {
                Title = "Нормальное распределение",
                Values = new ChartValues<ObservablePoint>(),
                MaxColumnWidth = 60
            });
            SeriesCollectionList.Add(new LineSeries
            {
                Title = "Нормальное распределение",
                Values = new ChartValues<ObservablePoint>()
            });
        }

        private double GenerateNormalDistributionRandomValue(double mean, double stDev)
        {
            var u1 = 1.0 - _rand.NextDouble();
            var u2 = 1.0 - _rand.NextDouble();
            var randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
            var randNormal = mean + stDev * randStdNormal;
            return randNormal;
        }

        private double GenerateUniformDistributionRandomValue(double mean, double stDev)
        {
            var a = mean - Math.Sqrt(3) * Math.Sqrt(stDev);
            var b = mean + Math.Sqrt(3) * Math.Sqrt(stDev);
            return a + _rand.NextDouble() * (b - a);
        }

        private double GenerateExponentialDistributionRandomValue(double mean, double? stDev)
        {
            var lambda = 1 / mean;
            var u = _rand.NextDouble();
            return (-lambda) * Math.Log(1-u);
        }

        //Нормальное распределение
        private void CreateNormalDistributionPlot(int samplingCount, double mean, double stDev, int beansCount)
        {
            _frequencies = new List<int>(new int[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++) //выборка по нормальному
            {
                _randoms[i] = GenerateNormalDistributionRandomValue(mean, stDev);
            }

            _randoms.Sort();
            var minVal = _randoms[0];
            var maxVal = _randoms[^1];
            var incr = (Math.Abs(minVal) + Math.Abs(maxVal)) / beansCount;

            var t = 0;

            var vals = new List<double>();
            double totalOfHeightBars = 0;
            for (var i = minVal; i < maxVal - incr; i += incr)
            {
                var count = _randoms.Count(c => c >= i && c <= i + incr);
                _frequencies[t] = count;
                totalOfHeightBars += _frequencies[t];

                SeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[t]));
                vals.Add(i);
                t++;
            }

            vals.ForEach(val =>
            {
                var yMean = (2 * val + incr) / 2;
                var y = NormalDistributionDensityFunction(stDev, mean, yMean);
                SeriesCollectionList?[1].Values.Add(new ObservablePoint(val, y * incr * totalOfHeightBars));
            });
        }

        //Равномерное распределение
        private void CreateUniformDistributionPlot(int samplingCount, double mean, double stDev, int beansCount)
        {
            _frequencies = new List<int>(new int[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++) //выборка по нормальному
            {
                _randoms[i] = GenerateUniformDistributionRandomValue(mean, stDev);
            }

            _randoms.Sort();
            var minVal = _randoms[0];
            var maxVal = _randoms[^1];
            var incr = (Math.Abs(minVal) + Math.Abs(maxVal)) / beansCount;

            var t = 0;

            var vals = new List<double>();
            double totalOfHeightBars = 0;
            for (var i = minVal; i < maxVal - incr; i += incr)
            {
                var count = _randoms.Count(c => c >= i && c <= i + incr);
                _frequencies[t] = count;
                totalOfHeightBars += _frequencies[t];

                SeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[t]));
                vals.Add(i);
                t++;
            }

            vals.ForEach(v => SeriesCollectionList?[1].Values.Add(new ObservablePoint(v, incr * totalOfHeightBars / (maxVal - minVal))));

        }

        private void CreateExponentialDistributionPlot(int samplingCount, double mean, double stDev, int beansCount)
        {
            _frequencies = new List<int>(new int[beansCount]);
            _randoms = new List<double>(new double[samplingCount]);

            for (var i = 0; i < samplingCount; i++) //выборка по нормальному
            {
                _randoms[i] = GenerateExponentialDistributionRandomValue(mean, stDev);
            }

            _randoms.Sort();
            var minVal = _randoms[0];
            var maxVal = _randoms[^1];
            var incr = (Math.Abs(minVal) + Math.Abs(maxVal)) / beansCount;

            var t = 0;

            var vals = new List<double>();
            double totalOfHeightBars = 0;
            for (var i = minVal; i < maxVal - incr; i += incr)
            {
                var count = _randoms.Count(c => c >= i && c <= i + incr);
                _frequencies[t] = count;
                totalOfHeightBars += _frequencies[t];

                SeriesCollectionList?[0].Values.Add(new ObservablePoint(i, _frequencies[t]));
                vals.Add(i);
                t++;
            }

            
            vals.ForEach(v =>
            {
                var yMean = (2 * v + incr) / 2;
                var y = ExponentialDistributionDensityFunction(1 / mean, yMean);
                SeriesCollectionList?[1].Values.Add(new ObservablePoint(v,
                    incr * totalOfHeightBars * y));
            });

        }
        private static double NormalDistributionDensityFunction(double sigm, double mean, double x)
        {
            return (1 / (sigm * Math.Sqrt(2 * Math.PI))) * Math.Pow(Math.E, (-1 * ((Math.Pow(x - mean, 2)) / (2 * sigm * sigm))));
        }

        private static double ExponentialDistributionDensityFunction(double lambda, double x)
        {
            return lambda * Math.Exp(-lambda * x);
        }

    }
}
