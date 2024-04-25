namespace DataStructs;

public class ArrayOfBorders
{
    public List<List<int>> Arr;

    public int Length { get; set; }

    public List<int> this[int i] => Arr[i];

    public MyEnumerator GetEnumerator() => new(this);

    public class MyEnumerator
    {  
        int nIndex;  
        ArrayOfBorders collection;  
        public MyEnumerator(ArrayOfBorders coll)
        {  
            collection = coll;  
            nIndex = -1;  
        }  
    
        public bool MoveNext()
        {  
            nIndex++;  
            return nIndex < collection.Arr.Count;  
        }  
    
        public List<int> Current => collection.Arr[nIndex];
    }

    private void SortList()
    {
        for (int i = 1; i < Length; i++)
            for (int j = 0; j < Length - i; j++)
                if (Arr[j][0] < Arr[j + 1][0])
                    (Arr[j], Arr[j + 1]) = (Arr[j + 1], Arr[j]);
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


    public void Add(List<int> borders) => Arr.Add(borders);

    public ArrayOfBorders(object o) => Arr = [];
}
