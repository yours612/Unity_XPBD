using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class Triangle
{
	public int Vertex1 { get; }
	public int Vertex2 { get; }
	public int Vertex3 { get; }

	public Triangle(int v1, int v2, int v3)
	{
		Vertex1 = v1;
		Vertex2 = v2;
		Vertex3 = v3;
	}
}
public class GetInformation : MonoBehaviour
{

	int[] edgesId;
	float[] verticesId;
	int[] oppositeEdgesIs;
	Mesh mymesh;
    void Start()
    {
		string filePath = "C:\\Users\\GMJ\\Desktop\\oppo.txt";
		StreamWriter sw;
		FileInfo fi = new FileInfo(filePath);
		sw = fi.AppendText();

		mymesh = new Mesh();

        int[] triangles = mymesh.tetSurfaceTriIds;
		Vector3[] vertices = mymesh.GetVPos();

		//��ȡ���б�
		edgesId = GetEdgeId(triangles);

		oppositeEdgesIs = GetOppositeEdges(triangles, edgesId);
		Debug.Log(edgesId.Length);
		Debug.Log(oppositeEdgesIs.Length);
		for (int i = 0; i < oppositeEdgesIs.Length; ++i) {
			//sw.Write(oppositeEdgesIs[i].ToString() + ",");
		}
		//sw.Close();
		//sw.Dispose();

		//��ȡ�����
		//verticesId = GetVerticesId(vertices);
		//Debug.Log(verticesId);
	}

	public int[] GetEdgeId(int[] triangles)
	{
		Dictionary<string, int> edges = new Dictionary<string, int>();

		for (int i = 0; i < triangles.Length; i += 3)
		{
			int vertex1 = triangles[i];
			int vertex2 = triangles[i + 1];
			int vertex3 = triangles[i + 2];

			// ��ӱ�1
			AddEdge(edges, vertex1, vertex2);
			// ��ӱ�2
			AddEdge(edges, vertex2, vertex3);
			// ��ӱ�3
			AddEdge(edges, vertex3, vertex1);
		}
		// �� Dictionary ת��Ϊ����
		int[] edgeArray = new int[edges.Count];
		edges.Values.CopyTo(edgeArray, 0);
		return edgeArray;
	}
	void AddEdge(Dictionary<string, int> edgeDict, int vertex1, int vertex2)
	{
		// ����ߵ��������������Ƿ���
		int smaller = Mathf.Min(vertex1, vertex2);
		int larger = Mathf.Max(vertex1, vertex2);
		string edgeKey = smaller.ToString() + "_" + larger.ToString();

		// ����������Ѿ��� Dictionary �У��Ͳ������
		if (!edgeDict.ContainsKey(edgeKey))
		{
			edgeDict.Add(edgeKey, smaller);
			edgeDict.Add(edgeKey + ".", larger);
		}
	}

	public float[] GetVerticesId(Vector3[] vertices) {
		float[] pos;
		pos = new float[vertices.Length * 3];

		for (int i = 0, j = 0; i < vertices.Length; ++i) {
			pos[j++] = vertices[i].x;
			pos[j++] = vertices[i].y;
			pos[j++] = vertices[i].z;
		}
		return pos;
	}

	public int[] GetOppositeEdges(int[] triangle, int[] edgesId) {
		List<int> oE = new List<int>(); 

		for (int j = 0; j < edgesId.Length; j += 2) {
			List<int> oV= new List<int>();
			for (int i = 0; i < triangle.Length; i += 3)
			{
				Triangle triangle1 = new Triangle(triangle[i], triangle[i + 1], triangle[i + 2]);
				int x = IsEdgeInTriangle(edgesId[j], edgesId[j + 1], triangle1);
				if (x != -1)
				{
					oV.Add(x);
				}
			}
			if (oV.Count != 0) {
				oE.Add(oV[0]);
				oE.Add(oV[1]);
			}
		}
		return oE.ToArray();
	}

	public int IsEdgeInTriangle(int edge_v1, int edge_v2, Triangle triangle)
	{
		if ((edge_v1 == triangle.Vertex1 || edge_v1 == triangle.Vertex2 || edge_v1 == triangle.Vertex3)
			&& (edge_v2 == triangle.Vertex1 || edge_v2 == triangle.Vertex2 || edge_v2 == triangle.Vertex3))
		{

			if (edge_v1 != triangle.Vertex1 && edge_v2 != triangle.Vertex1)
			{
				return triangle.Vertex1;
			}
			else if (edge_v1 != triangle.Vertex2 && edge_v2 != triangle.Vertex2)
			{
				return triangle.Vertex2;
			}
			else
			{
				return triangle.Vertex3;
			}
		}
		else {
			return -1;
		}
		
	}

}
