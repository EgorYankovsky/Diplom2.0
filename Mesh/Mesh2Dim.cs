namespace Grid;
using System.Collections.Immutable;
using DataStructs;

public class Mesh2Dim : Mesh
{
    public override int NodesAmountTotal 
    { 
        get => NodesAmountR * NodesAmountZ;
    }

    public int NodesAmountR
    { 
        get => nodesR.Count; 
    }

    internal List<int> nodesR_Refs;

    internal ImmutableArray<double> NodesRWithoutFragmentation { get; set; }

    internal string? infoAboutR;

    public int NodesAmountZ 
    { 
        get => nodesZ.Count;
    }

    internal List<int> nodesZRefs;

    public ImmutableArray<double> NodesZWithoutFragmentation { get; set; }

    internal string? infoAboutZ;
    

    public Mesh2Dim()
    {
        borders = new();
        Elems = new();
        nodesZ = new();
        nodesR = new();
        nodesR_Refs = new();
        nodesZRefs = new();
        mu0 = new();
        sigma = new();
    }
}
