using System.Collections;

namespace DataStructs;

public class ArrayOfBorders : IEnumerable
{
    public List<List<int>> Arr;

    public int Length { get; set; }

    public List<int> this[int i] => Arr[i];

    public IEnumerator<List<int>> GetEnumerator() => Arr.GetEnumerator();

    private void SortList()
    {
        for (int i = 1; i < Length; i++)
            for (int j = 0; j < Length - i; j++)
                if (Arr[j][0] < Arr[j + 1][0])
                    (Arr[j], Arr[j + 1]) = (Arr[j + 1], Arr[j]);
    }

    public void Add(List<int> arr) => Arr.Add(arr);

    IEnumerator IEnumerable.GetEnumerator()
    {
        throw new NotImplementedException();
    }

    public ArrayOfBorders(string path)
    {
        Arr = [];
        var data = File.ReadAllText(path).Split("\n");
        Length = int.Parse(data[0]);
        for (int i = 0; i < Length; i++)
            Arr.Add([.. data[i + 1].Split(" ").Select(int.Parse)]);
        SortList();
    }

    public ArrayOfBorders(string path, int f)
    {
        Arr = [];
        var data = File.ReadAllText(path).Split("\n");
        Length = int.Parse(data[0]);
        for (int i = 0; i < Length; i++)
            Arr.Add([.. data[i + 1].Split(" ").Select(int.Parse)]);
    }
}
