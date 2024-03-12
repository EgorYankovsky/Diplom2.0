using System.Diagnostics;

namespace Processor;

public static class Postprocessor
{
    private static readonly string _arg1 = Path.GetFullPath("../../../../Data/Subtotals/Points.dat");

    private static readonly string _arg2 = Path.GetFullPath("../../../../Data/Output/");

    private static readonly string _arg3 = Path.GetFullPath("../../../../Drawer/Pictures/");

    private static readonly string _source = "..\\..\\..\\..\\Drawer\\source.py";

    private static readonly string _fileName = "python";

    private static Process? process;

    public static int ShowAnswer()
    {
        process = new();
        process.StartInfo.Arguments = $"{_source} {_arg1} {_arg2} {_arg3}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }
}
