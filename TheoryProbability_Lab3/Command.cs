using System;
using System.Windows.Input;

namespace TheoryProbability_Lab3
{
    public class Command : ICommand
    {
        private readonly Action<object> _execute;
        private Func<object, bool> _canExecute;

        public Command(Action<object> execute, Func<object, bool> canExecute)
        {
            _execute = execute;
            _canExecute = canExecute;
        }
        public event EventHandler CanExecuteChanged;

        public bool CanExecute(object parameter)
        {
            return true;
        }

        public void Execute(object parameter)
        {
            _execute(parameter);
        }
    }
}
