namespace Manager;

public static class FolderManager
{
    public static void ClearFolder(string path)
    {
        DirectoryInfo di = new(path);
        foreach (FileInfo file in di.GetFiles())
            file.Delete(); 
    }
}
