using System.Collections;

namespace DataStructs;

public class ArrayOfElems : IEnumerable
{
    public List<Elem> elems;

    public int Length { get => elems.Count; set {} }

    public IEnumerator GetEnumerator() => elems.GetEnumerator();

    public Elem this[int i] => elems[i];

    public List<int> GetLinks(int i) => elems[i].Arr;

    public double GetMu(int i) => elems[i].mu;

    public double GetSigma(int i) => elems[i].sigma;

    public ArrayOfElems(string path)
    {
        elems = [];
        var data = File.ReadAllText(path).Split("\n");
        int length = int.Parse(data[0]);
        for (int i = 0; i < length; i++)
        {
            var info = data[1 + i].Split(" ");
            elems.Add(new Elem(0, int.Parse(info[0]), int.Parse(info[1]), int.Parse(info[2]),
                               int.Parse(info[3]), double.Parse(info[4]), double.Parse(info[5])));
        }
    }

    public ArrayOfElems(string path, int num)
    {
        elems = [];
        var data = File.ReadAllText(path).Split("\n");
        int length = int.Parse(data[0]);
        for (int i = 0; i < length; i++)
        {
            var info = data[1 + i].Split(" ");
            List<int> linksArr = info[1 .. 13].Select(int.Parse).ToList();
            elems.Add(new Elem(int.Parse(info[0]), linksArr, double.Parse(info[13]), double.Parse(info[14])));
        }
    }

    public void Add(Elem elem) => elems.Add(elem);
}