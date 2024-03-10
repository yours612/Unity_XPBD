using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class _ReadFile : MonoBehaviour
{
    public List<int> tetId = new List<int>();
    public List<float> vertsPos = new List<float>();
    public List<int> tetEdgeIds = new List<int>();
    public List<int> tetSurfaceTriIds = new List<int>();

    public void ReadFiles()
    {
        ReadElefile();
        ReadNodeFile();
        ReadEdgeFile();
        ReadFaceFile();
    }

    List<string> SplitString(string str, char delimiter1)
    {
        List<string> tokens = new List<string>();
        StringReader reader = new StringReader(str);

        char delimiter2 = '\n';
        string line;

        while ((line = reader.ReadLine()) != null)
        {
            string[] subTokens = line.Split(delimiter1);
            foreach (string subToken in subTokens)
            {
                if (!string.IsNullOrEmpty(subToken))
                {
                    tokens.Add(subToken);
                }
            }
        }

        return tokens;
    }
    void ReadElefile()
    {
        string filePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/HighDannang/High.1.ele";
        //string filePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/lowDannang/Low.1.ele";//低模的
        //List<int> tetId = new List<int>();

        using (StreamReader file = new StreamReader(filePath))
        {
            if (file == null)
            {
                return;
            }

            string fileContent = file.ReadToEnd();
            file.Close();

            List<string> strings = SplitString(fileContent, ' ');
            int tetNumber = int.Parse(strings[0]);

            tetId.Capacity = tetNumber * 4;

            for (int tet = 0; tet < tetNumber; ++tet)
            {
                tetId.Add(int.Parse(strings[tet * 5 + 4]));
                tetId.Add(int.Parse(strings[tet * 5 + 5]));
                tetId.Add(int.Parse(strings[tet * 5 + 6]));
                tetId.Add(int.Parse(strings[tet * 5 + 7]));
            }
        }

        //Debug.Log(tetId.Count);
    }
    void ReadNodeFile()
    {
        string filePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/HighDannang/High.1.node";
        //string filePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/lowDannang/Low.1.node";//低模的

        using (StreamReader file = new StreamReader(filePath))
        {
            if (file == null)
            {
                return;
            }

            string fileContent = file.ReadToEnd();
            file.Close();

            List<string> strings = SplitString(fileContent, ' ');
            int xNumber = int.Parse(strings[0]);

            vertsPos.Capacity = xNumber * 3;

            for (int i = 0; i < xNumber; ++i)
            {
                vertsPos.Add(float.Parse(strings[i * 4 + 5]));
                vertsPos.Add(float.Parse(strings[i * 4 + 6]));
                vertsPos.Add(float.Parse(strings[i * 4 + 7]));
            }
        }

        //Debug.Log(vertsPos.Count);
    }
    void ReadEdgeFile()
    {
        string edgeFilePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/HighDannang/High.1.edge";
        //string edgeFilePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/lowDannang/Low.1.edge";//低模的

        using (StreamReader file = new StreamReader(edgeFilePath))
        {
            if (file == null)
            {
                return;
            }

            string fileContent = file.ReadToEnd();
            file.Close();

            List<string> strings = SplitString(fileContent, ' ');

            int eNum = int.Parse(strings[0]);
            tetEdgeIds.Capacity = eNum * 2;

            for (int i = 0; i < eNum; ++i)
            {
                tetEdgeIds.Add(int.Parse(strings[i * 4 + 3]));
                tetEdgeIds.Add(int.Parse(strings[i * 4 + 4]));
            }
        }
        //Debug.Log(tetEdgeIds.Count);
    }
    void ReadFaceFile()
    {
        string faceFilePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/HighDannang/High.1.face";
        //string faceFilePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/lowDannang/Low.1.face";//低模的

        using (StreamReader file = new StreamReader(faceFilePath))
        {
            if (file == null)
            {
                return;
            }

            string fileContent = file.ReadToEnd();
            file.Close();

            List<string> strings = SplitString(fileContent, ' ');

            int fNum = int.Parse(strings[0]);
            tetSurfaceTriIds.Capacity = fNum * 3;

            for (int i = 0; i < fNum; ++i)
            {
                tetSurfaceTriIds.Add(int.Parse(strings[i * 5 + 3]));
                tetSurfaceTriIds.Add(int.Parse(strings[i * 5 + 4]));
                tetSurfaceTriIds.Add(int.Parse(strings[i * 5 + 5]));
            }
        }
        //Debug.Log(tetSurfaceTriIds.Count);
    }

}

