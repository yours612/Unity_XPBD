using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class MSM : MonoBehaviour
{
    UnityEngine.Mesh graphMesh;
    MeshCollider meshCollider;

    List<float> vertsPos = new List<float>();
    float[] vertspos = { };
    List<int> tetSurfaceTriIds = new List<int>();

    Vector3[] Xpos;
    public Vector3[] originXpos;
    int[] face;
    int[] edgeId;
    List<int[]> edges = new List<int[]>();
    float[] edgeLengths;

    public float vk;
    public float ks;
    public float kd;
    public float gravity;
    public float dt;

    Vector3[] vel;
    Vector3[] force;

    public int pickVertex = -1;
    Vector3 pickVertexPos = new Vector3(0, 0, 0);

    void Start()
    {
        ReadNodeFile();
        ReadFaceFile();
        vertspos = vertsPos.ToArray();
        face = tetSurfaceTriIds.ToArray();

        graphMesh = GetComponent<MeshFilter>().mesh;
        meshCollider = GetComponent<MeshCollider>();

        originXpos = new Vector3[vertspos.Length/3];
        Xpos = new Vector3[vertspos.Length / 3];
        for (int i = 0, j = 0; i < vertspos.Length; i += 3)
        {
            Vector3 vertex;
            vertex.x = vertspos[i];
            vertex.y = vertspos[i + 1] ;
            vertex.z = vertspos[i + 2] ;
            Xpos[j] = vertex;
            originXpos[j++] = vertex;
        }

        GetEdge();
        edgeId = new int[edges.Count * 2];
        for (int i = 0, j = 0; i < edges.Count; ++i) {
            edgeId[j++] = edges[i][0];
            edgeId[j++] = edges[i][1];
        }

        
        edgeLengths = new float[edgeId.Length / 2];
        for (int i = 0; i < edgeLengths.Length; i++)
        {
            int id0 = edgeId[2 * i];
            int id1 = edgeId[2 * i + 1];
            edgeLengths[i] = Vector3.Distance(Xpos[id0], Xpos[id1]);
        }
        
        vel = new Vector3[Xpos.Length];
        force = new Vector3[Xpos.Length];
        for (int i = 0; i < Xpos.Length; i++) {
            vel[i] = new Vector3(0.0f, 0, 0);
            force[i] = new Vector3(0.0f, 0, 0);
        }

        graphMesh.vertices = Xpos;
        graphMesh.triangles = face;
        graphMesh.RecalculateNormals();

       //GetComponent<YanCheZuo.TestMeshViewer>().enabled = true;
    }

    public void initialize() {
        ReadNodeFile();
        ReadFaceFile();
        vertspos = vertsPos.ToArray();
        face = tetSurfaceTriIds.ToArray();

        originXpos = new Vector3[vertspos.Length / 3];
        Xpos = new Vector3[vertspos.Length / 3];
        for (int i = 0, j = 0; i < vertspos.Length; i += 3)
        {
            Vector3 vertex;
            vertex.x = vertspos[i];
            vertex.y = vertspos[i + 1];
            vertex.z = vertspos[i + 2];
            Xpos[j] = vertex;
            originXpos[j++] = vertex;
        }
    }

    void Update()
    {
        for (int i = 0; i < Xpos.Length; ++i)
        {
            force[i].y += gravity;
            vel[i] *= kd;            
        }

        StartCoroutine(OnMouseDown());
        if (Input.GetMouseButtonDown(0))
        {
            PickVertex();
        }
        if (pickVertex != -1)
        {
            Xpos[pickVertex] = pickVertexPos;
        }

        //自身点约束
        for (int i = 0; i < Xpos.Length; ++i)
        {
            Vector3 director = Xpos[i] - originXpos[i];
            float len = director.magnitude;
            director = director / len;

            if (Mathf.Abs(len) > 0.01f)
            {
                float detaLen = len;
                if (Mathf.Abs(detaLen) > 3.0f)
                {
                    Xpos[i] = originXpos[i] + director * 3;
                    detaLen = 3;
                }
                force[i] += -vk * director * detaLen;
            }
        }
        //边约束
        for (int i = 0; i < edgeLengths.Length; i++) 
        {
            int id0 = edgeId[2 * i];
            int id1 = edgeId[2 * i + 1];

            Vector3 director = Xpos[id0] - Xpos[id1];
            float len = director.magnitude;
            director = director / len;

            float detaLen = len - edgeLengths[i];
            if (Mathf.Abs(detaLen) > 0.01f)
            {
                force[id0] += -ks * director * detaLen;
                force[id1] -= -ks * director * detaLen;
            }
        }
        

        for (int i = 0; i < Xpos.Length; ++i)
        {
            vel[i] += force[i] * dt;
            Xpos[i] += vel[i] * dt;
            force[i] = new Vector3(0, 0, 0);
        }

        
        graphMesh.vertices = Xpos;
        graphMesh.triangles = face;
        graphMesh.RecalculateNormals();
        meshCollider.sharedMesh = graphMesh;
    }
    

    void GetEdge() {
        for (int i = 0; i < face.Length; i += 3) {
            int v1 = face[i];
            int v2 = face[i + 1];
            int v3 = face[i + 2];
            AddEdge(v1, v2);
            AddEdge(v2, v3);
            AddEdge(v1, v3);
        }
    }
    void AddEdge(int index1, int index2) {
        int[] newEdge = { index1, index2 }; 
        if (edges.Find(t => (t[0] == index1 && t[1] == index2)
                     || (t[0] == index2 && t[1] == index1)) == null)
        {
                edges.Add(newEdge);
        }
        
    }


    IEnumerator OnMouseDown()
    {
        if (pickVertex != -1)
        {
            Vector3 targetScreenPos = Camera.main.WorldToScreenPoint(Xpos[pickVertex]);
            //var offset = Xpos[pickVertex] - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, targetScreenPos.z));

            while (Input.GetMouseButton(0))
            {
                Vector3 mousePos = new Vector3(Input.mousePosition.x, Input.mousePosition.y, targetScreenPos.z);
                var targetPos = Camera.main.ScreenToWorldPoint(mousePos);// + offset;
                pickVertexPos = targetPos;
                yield return new WaitForFixedUpdate();//循环执行
            }
        }
        pickVertex = -1;
    }
    void PickVertex()
    {
        RaycastHit hit;
        Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

        if (Physics.Raycast(ray, out hit))
        {
            if (hit.collider == null)
                return;
            MeshCollider collider = (MeshCollider)hit.collider;
            if (collider.gameObject.CompareTag("liver"))
            {
                //Debug.Log(hit.collider.gameObject.name);
                UnityEngine.Mesh mesh0 = collider.sharedMesh;
                Vector3[] vertices = mesh0.vertices;
                int[] triangles = mesh0.triangles;
                Vector3 p0 = hit.transform.TransformPoint(vertices[triangles[hit.triangleIndex * 3]]);
                Vector3 p1 = hit.transform.TransformPoint(vertices[triangles[hit.triangleIndex * 3 + 1]]);
                Vector3 p2 = hit.transform.TransformPoint(vertices[triangles[hit.triangleIndex * 3 + 2]]);
                pickVertex = MinDistanceMouse(p0, p1, p2, triangles, hit.triangleIndex);
            }
        }
    }
    int MinDistanceMouse(Vector3 v1, Vector3 v2, Vector3 v3, int[] triangles, int index)
    {
        Vector3 mouseposition = Input.mousePosition;
        float distance1 = Vector3.Distance(v1, mouseposition);
        float distance2 = Vector3.Distance(v2, mouseposition);
        float distance3 = Vector3.Distance(v3, mouseposition);
        if (distance1 < distance2 && distance1 < distance3)
            return triangles[index * 3];
        if (distance2 < distance1 && distance2 < distance3)
            return triangles[index * 3 + 1];
        if (distance3 < distance1 && distance3 < distance2)
            return triangles[index * 3 + 2];
        return -1;
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
    void ReadNodeFile() {
        string filePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/LiverPlane/liverPlane_vert.txt";
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
                vertsPos.Add(float.Parse(strings[i * 4 + 2]));
                vertsPos.Add(float.Parse(strings[i * 4 + 4]) * -1);
                vertsPos.Add(float.Parse(strings[i * 4 + 3]));
            }
        }

    }
    void ReadFaceFile()
    {
        string faceFilePath = "E:/tetgen/tetgen-master/vs/tetgen/Debug/LiverPlane/liverPlane_face.txt";
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
                List<string> strings1 = SplitString(strings[i * 4 + 2], '/');
                tetSurfaceTriIds.Add(int.Parse(strings1[0]) - 1);

                strings1 = SplitString(strings[i * 4 + 3], '/');
                tetSurfaceTriIds.Add(int.Parse(strings1[0]) - 1);

                strings1 = SplitString(strings[i * 4 + 4], '/');
                tetSurfaceTriIds.Add(int.Parse(strings1[0]) - 1);
            }
        }
        //Debug.Log(tetSurfaceTriIds.Count);
    }
}
