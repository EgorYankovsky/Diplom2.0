using System.Collections;

namespace Manager;

public static class FolderManager
{
    public static void ClearFolder(string path)
    {
        DirectoryInfo di = new(path);
        foreach (FileInfo file in di.GetFiles())
            file.Delete();
    }

    public static void ClearFolders(IEnumerable paths)
    {
        foreach (string path in paths)
            ClearFolder(path);
    }
}
