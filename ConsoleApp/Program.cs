using Project;
using System.Globalization;
using Dumpify;
using System.Diagnostics;
using Solver;
using Manager;
using Processor;
using Solution;
using MathObjects;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
DumpConfig.Default.ColorConfig.ColumnNameColor = DumpColor.FromHexString("#FF0000");
DumpConfig.Default.ColorConfig.PropertyValueColor = DumpColor.FromHexString("#00FF00");
DumpConfig.Default.ColorConfig.NullValueColor = DumpColor.FromHexString("#0000FF");

string CalculationArea = Path.GetFullPath("../../../../Data/Input/Info.dat");
string BordersInfo = Path.GetFullPath("../../../../Data/Input/Borders.dat");
string AnswerPath = Path.GetFullPath("../../../../Data/Output/");
string SubtotalsPath = Path.GetFullPath("../../../../Data/Subtotals/");
string TimePath = Path.GetFullPath("../../../../Data/Input/Time.dat");
string PicturesPath = Path.GetFullPath("../../../../Drawer/Pictures/");

// Pre-processor.
FolderManager.ClearFolder(AnswerPath);
FolderManager.ClearFolder(PicturesPath);
FolderManager.ClearFolder(SubtotalsPath);

FEM2D myFEM2D = new();
myFEM2D.ReadData(CalculationArea, BordersInfo, TimePath);   
myFEM2D.ConstructMesh();
myFEM2D.SubmitGeneratedData();
myFEM2D.SetSolver(new LOS());
myFEM2D.Solve();
myFEM2D.WriteData(AnswerPath);
myFEM2D.GenerateVectorEphi();

FEM3D myFEM3D = new(myFEM2D);
myFEM3D.ConstructMesh(myFEM2D);
myFEM3D.GenerateExy(myFEM2D);



int a = Postprocessor.ShowAnswer();
Debug.WriteLine($"Postprocessor finished with code: {a}");