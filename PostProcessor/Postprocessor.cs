using System.Diagnostics;

namespace Processor;

public static class Postprocessor
{
    private static readonly string _arg1 = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Points.poly");

    private static readonly string _arg2 = Path.GetFullPath("../../../../Data/Output/A_phi/Answer/");

    private static readonly string _arg3 = Path.GetFullPath("../../../../Drawer/Pictures/A_phi/");

    private static readonly string _arg4 = Path.GetFullPath("../../../../Data/Output/E_phi/Answer/");

    private static readonly string _arg5 = Path.GetFullPath("../../../../Drawer/Pictures/E_phi/");

    private static readonly string _source = "..\\..\\..\\..\\Drawer\\source.py";

    private static readonly string _sourceAnimation = "..\\..\\..\\..\\Drawer\\animation_creator.py";


    private static readonly string _fileName = "python";

    private static Process? process;

    public static int DrawA_phi()
    {
        process = new();
        process.StartInfo.Arguments = $"{_source} {_arg1} {_arg2} {_arg3}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        int animStatus = DrawA_phi_gif();
        return process.ExitCode & animStatus;
    }

    public static int DrawE_phi()
    {
        process = new();
        process.StartInfo.Arguments = $"{_source} {_arg1} {_arg4} {_arg5}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        int animStatus = DrawE_phi_gif();
        return process.ExitCode & animStatus;
    }

    private static int DrawA_phi_gif()
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceAnimation} {_arg1} {_arg2} {_arg3}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }

    private static int DrawE_phi_gif()
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceAnimation} {_arg1} {_arg4} {_arg5}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }
}
