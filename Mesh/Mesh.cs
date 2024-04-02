namespace Grid;

public abstract class Mesh
{
    public abstract int ElemsAmount { get; set; }
 
    internal int bordersAmount;
 
    internal List<List<int>>? borders;

    public List<double>? mu0;

    public List<double>? sigma;

    public List<double>? nodesR;
    
    public List<double>? nodesX;
    
    public List<double>? nodesY;
    
    public List<double>? nodesZ;

    public List<List<int>>? Elems;

    public abstract int NodesAmountTotal { get; }

    public int NodesAmountX => nodesX.Count;
    public int NodesAmountY => nodesY.Count;
    public int NodesAmountZ => nodesZ.Count;
}