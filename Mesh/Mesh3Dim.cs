namespace Grid;
using System.Collections.Immutable;
using DataStructs;

public class Mesh3Dim : Mesh
{
    public override int NodesAmountTotal 
    { 
        get => NodesAmountX * NodesAmountY * NodesAmountZ;
    }

    public override int ElemsAmount
    {
        get => (NodesAmountX - 1) * (NodesAmountY - 1) * (NodesAmountZ - 1);
        set => ElemsAmount = value;
    }

    internal List<int> nodesXRefs;

    internal ImmutableArray<double> NodesXWithoutFragmentation { get; set; }

    internal string? infoAboutX;

    internal List<int> nodesYRefs;

    internal ImmutableArray<double> NodesYWithoutFragmentation { get; set; }

    internal string? infoAboutY;

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