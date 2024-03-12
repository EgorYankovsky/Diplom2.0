namespace Grid;
using System.Collections.Immutable;
using DataStructs;

public class Mesh3Dim : Mesh
{
    public override int NodesAmountTotal 
    { 
        get => NodesAmountX * NodesAmountY * NodesAmountZ;
    }

    public int NodesAmountX
    { 
        get => nodesX.Count; 
    }

    internal List<int> nodesXRefs;

    internal ImmutableArray<double> NodesXWithoutFragmentation { get; set; }

    internal string? infoAboutX;


    public int NodesAmountY
    { 
        get => nodesY.Count; 
    }

    internal List<int> nodesYRefs;

    internal ImmutableArray<double> NodesYWithoutFragmentation { get; set; }

    internal string? infoAboutY;


    public int NodesAmountZ 
    { 
        get => nodesZ.Count;
    }

    internal List<int> nodesZRefs;

    public ImmutableArray<double> NodesZWithoutFragmentation { get; set; }

    internal string? infoAboutZ;
    
    public Mesh3Dim()
    {
        borders = new();
        Elems = new();
        nodesZ = new();
        nodesX = new();
        nodesY = new();
        nodesXRefs = new();
        nodesZRefs = new();
        nodesYRefs = new();
        mu0 = new();
        sigma = new();
    }
}