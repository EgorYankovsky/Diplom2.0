using DataStructs;

namespace Grid;

public abstract class Mesh(List<Elem> elems)
{
    public abstract int ElemsAmount { get; set; }
 
    internal int bordersAmount;

    public List<Elem> Elems = elems;

    public abstract int NodesAmountTotal { get; }
}