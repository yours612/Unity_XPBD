using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public struct _Matrix3x3
{
	public float M11, M12, M13;
	public float M21, M22, M23;
	public float M31, M32, M33;

	public static float Determinant(_Matrix3x3 matrix)
	{
		return matrix.M11 * matrix.M22 * matrix.M33 +
			   matrix.M12 * matrix.M23 * matrix.M31 +
			   matrix.M13 * matrix.M21 * matrix.M32 -
			   matrix.M31 * matrix.M22 * matrix.M13 -
			   matrix.M32 * matrix.M23 * matrix.M11 -
			   matrix.M33 * matrix.M21 * matrix.M12;
	}
}

public class Low2High_XPBD : MonoBehaviour
{
	public int test = 0;

	public float K = 0.001f;
	public bool openEdgeConstrian = true;
	public float edgeK = 0.001f;
	public bool openVolumeConstrian = false;
	public float volumeK = 1.0f;
	public float edgeCompliance;
	public float volCompliance;
	public float gravity;
	public float damping = 0.99f;


	//-----------低模：四面体tet、顶点vert、边表edge      高模：顶点vert， 三角面tri-----------------

	//低模对应高模的顶点的连线向量----高模的顶点与最近的低模点相连
	Dictionary<int, Dictionary<int, Vector3>> Map = new Dictionary<int, Dictionary<int, Vector3>>();
	//Dictionary<int, Dictionary<int[], Vector3>> _Map = new Dictionary<int, Dictionary<int[], Vector3>>();
	Dictionary<int, Dictionary<int, Vector3>> Map_Second = new Dictionary<int, Dictionary<int, Vector3>>();
	Dictionary<int, int> mapForPick = new Dictionary<int, int>();

	//低模
	int Low_numVerts = 28;
	int Low_numTets;
	float[] Low_X = { -0.88251f, -2.545111f, 1.715908f, 1.461392f, -0.667046f, -1.083844f, 
		1.906528f, 1.307936f, -1.708346f, -0.638421f, 2.324212f, -0.506435f, 2.524736f, 0.32408f, -2.570727f, 
		0.286067f, -0.115257f, 0.688311f, -2.280253f, -1.168516f, 4.136124f, 1.508731f, 0.474631f, -0.094702f, 
		-2.95919f, -2.485795f, 1.446286f, 0.213531f, 2.019937f, 0.180592f, -0.962037f, 1.59776f, -1.260553f,
		0.874915f, 1.260822f, -1.820215f, 2.622734f, 0.629177f, -2.533818f, -3.291321f, 0.153379f, 0.432102f, 
		-3.721343f, -1.146807f, 2.777617f, -2.09509f, 1.370568f, 3.097165f, -2.99766f, 1.469091f, 2.283977f, 
		-1.329501f, 0.11014f, -1.278994f, 0.045444f, -0.389551f, -2.018918f, -1.91218f, -2.730891f, 3.13327f, 
		-0.492304f, -0.883444f, -0.830906f, 2.9442f, 0.325988f, -2.453042f, 3.071238f, 0.109163f, -1.774583f, 
		-0.596406f, -1.129216f, 3.592819f, -0.064425f, -0.908915f, -0.153355f, 2.995369f, 0.687531f, -2.433713f, 
		1.475488f, -0.20706f, -2.514508f, -0.339059f, 1.288017f, 0.973573f };
	int[] Low_face = {9, 2, 22, 19, 23, 0, 18, 20, 1, 0, 5, 24, 18, 11, 10, 13, 14, 8, 6, 19, 14, 22, 2, 25,
		27, 15, 3, 10, 2, 3, 7, 9, 22, 8, 24, 20, 3, 13, 10, 4, 26, 21, 12, 2, 11, 11, 26, 12, 9, 27, 3, 19, 0, 8,
		27, 5, 23, 23, 6, 15, 5, 7, 1, 17, 13, 8, 10, 17, 18, 5, 27, 7, 17, 10, 13, 3, 2, 9, 0, 23, 5, 10, 11, 2, 4,
		12, 26, 8, 14, 19, 15, 6, 16, 8, 0, 24, 3, 16, 13, 14, 16, 6, 3, 15, 16, 16, 14, 13, 17, 8, 20, 1, 20, 24, 26,
		18, 1, 26, 11, 18, 19, 6, 23, 1, 22, 26, 18, 17, 20, 21, 22, 25, 27, 23, 15, 24, 5, 1, 1, 7, 22, 12, 25, 2, 21,
		26, 22, 9, 7, 27, 4, 21, 12, 25, 12, 21};
	//用于旋转的三个向量组
	Vector3[] Low_OriginNormals;
	Vector3[] Low_CurNormals;
	//Vector3[] ForRotate_Line; // 判断扭转
	Dictionary<int, Dictionary<int, Vector3>> ForRotate_Line = new Dictionary<int, Dictionary<int, Vector3>>();// 判断扭转

	Vector3[] Low_Xpos;
	Vector3[] Low_prevX;
	Vector3[] vel;
	

	int[] Low_Tet;
	int[] Low_edgeIds;
	float[] restVol;
	float[] Low_edgeLengths;
	float[] invMass;

	//高模
	int[] High_Tris;
	int High_numVerts;
	float[] High_X;
	Vector3[] High_Xpos;
	float[] High_edgeLengths;
	int[] High_edgeIds;

	public Mesh mesh;

	int[,] volIdOrder = {
		{1, 3, 2},
		{0, 2, 3},
		{0, 3, 1},
		{0, 1, 2}
	};

	//图形显示
	UnityEngine.Mesh graphMesh;
	Bounds Bounds;

	MeshCollider meshCollider;
	Bounds bounds1;

	public int pickHigh_Vertex = -1;
	public int pickLow_Vertex = -1;
	Vector3 pickHigh_VertexPos = new Vector3(0, 0, 0);
	Vector3 pickLow_VertexPos = new Vector3(0, 0, 0);

	//fps
	public float fps;
	int frame;
	float duration;

	// 初始化参数
	void Start()
	{
		graphMesh = GetComponent<MeshFilter>().mesh;
		meshCollider = GetComponent<MeshCollider>();

		mesh = new Mesh();
		mesh.GetInfor();

		High_numVerts = mesh.vertsPos.Length / 3;
		High_X = mesh.vertsPos;

		Low_numTets = mesh.tetId.Length / 4;
		Low_Xpos = new Vector3[Low_numVerts];
		vel = new Vector3[Low_numVerts];
		Low_Tet = mesh.tetId;
		Low_edgeIds = mesh.tetEdgeIds;
		restVol = new float[Low_numTets];
		invMass = new float[Low_numVerts];
		Low_edgeLengths = new float[Low_edgeIds.Length / 2];

		Low_OriginNormals = new Vector3[Low_numVerts];
		Low_CurNormals = new Vector3[Low_numVerts];

		
		High_Xpos = new Vector3[High_numVerts];
		High_edgeIds = mesh.tetHighEdgeIds;
		High_edgeLengths = new float[High_edgeIds.Length / 2]; //高模再用距离约束--无效
		//edgeCompliance = 0;
		//volCompliance = 0;


		//初始化低模高模的顶点数组；
		for (int i = 0, j = 0; i < Low_X.Length; i += 3)
		{
			Vector3 vertex;
			vertex.x = Low_X[i] * 0.5f;
			vertex.y = Low_X[i + 1] * 0.5f;
			vertex.z = Low_X[i + 2] * 0.5f;
			Low_Xpos[j++] = vertex;
			//VirtualXpos[j++] = vertex;
		}
		Low_prevX = new Vector3[Low_Xpos.Length];
		for (int i = 0, j = 0; i < High_X.Length; i += 3) {
			Vector3 vertex;
			vertex.x = High_X[i] * 0.5f;
			vertex.y = High_X[i + 1] * 0.5f;
			vertex.z = High_X[i + 2] * 0.5f;
			High_Xpos[j++] = vertex;
		}
		//初始化低模顶点法线
		Low_OriginNormals = calculateLow_Normal();
		//初始化低模对应高模的顶点连线
		for (int i = 0; i < Low_numVerts; ++i)
		{ 
			Map[i] = new Dictionary<int, Vector3>();
			Map_Second[i] = new Dictionary<int, Vector3>();
			//_Map[i] = new Dictionary<int[], Vector3>();
			ForRotate_Line[i] = new Dictionary<int, Vector3>();
			int k;
			for (k = 0; k < Low_face.Length; k += 3)
			{
				if (Low_face[k] == i)
				{
					k = k + 1;
					break;
				}
				else if (Low_face[k+1] == i)
				{
					k = k;
					break;
				}
				else if (Low_face[k+2] == i)
				{
					k = k + 1;
					break;
				}
			}
			ForRotate_Line[i].Add(Low_face[k], Low_Xpos[i] - Low_Xpos[Low_face[k]]);
		}
		for (int j = 0; j < High_numVerts; ++j) {
			Vector3 pos = High_Xpos[j];
			float minDistance = float.MaxValue;
			float second_min = float.MaxValue;
			int closestLowPolyVertexIndex = -1;
			int second_closestV = -1;

			for (int i = 0; i < Low_numVerts; ++i) {
				float distance = Vector3.Distance(pos, Low_Xpos[i]);
				if (distance < minDistance)
				{
					second_min = minDistance;
					second_closestV = closestLowPolyVertexIndex;

					minDistance = distance;
					closestLowPolyVertexIndex = i;
				}
				else if (distance < second_min) {
					second_min = distance;
					second_closestV = i;
				}
			}

			if (closestLowPolyVertexIndex != -1) {
				mapForPick[j] = closestLowPolyVertexIndex;
				//Map[closestLowPolyVertexIndex].Add(j, Low_Xpos[closestLowPolyVertexIndex] - pos);
				//Map[closestLowPolyVertexIndex].Add(j, Vector3.Magnitude(Low_Xpos[closestLowPolyVertexIndex] - pos));
				Map[closestLowPolyVertexIndex].Add(j, Low_Xpos[closestLowPolyVertexIndex] - pos - Low_OriginNormals[closestLowPolyVertexIndex]);
				//_Map[closestLowPolyVertexIndex].Add(x, Low_Xpos[closestLowPolyVertexIndex] - pos - Low_OriginNormals[closestLowPolyVertexIndex]);
				Map_Second[second_closestV].Add(j, Low_Xpos[second_closestV] - pos - Low_OriginNormals[second_closestV]);
			}
		}
		
		
		

		High_Tris = new int[mesh.tetSurfaceTriIds.Length];
		for (int i = 0; i < mesh.tetSurfaceTriIds.Length; i += 3)
		{
			High_Tris[i] = mesh.tetSurfaceTriIds[i];
			High_Tris[i + 1] = mesh.tetSurfaceTriIds[i + 2];
			High_Tris[i + 2] = mesh.tetSurfaceTriIds[i + 1];
		}

		//所有数组置零
		for (int i = 0; i < Low_numVerts; ++i)
		{
			vel[i] = new Vector3(0.0f, 0.0f, 0.0f);
			restVol[i] = 0;
			invMass[i] = 0;
		}

		//初始化质量矩阵，剩余体积矩阵, 边长矩阵
		for (int i = 0; i < Low_numTets; i++)
		{
			float vol = getTetVolume(i);
			restVol[i] = vol;
			float pInvMass = vol > 0.0f ? 1.0f / (vol / 4.0f) : 0.0f;
			invMass[Low_Tet[4 * i]] += pInvMass;
			invMass[Low_Tet[4 * i + 1]] += pInvMass;
			invMass[Low_Tet[4 * i + 2]] += pInvMass;
			invMass[Low_Tet[4 * i + 3]] += pInvMass;
		}

		int[] value = { 1, 0, 8, 84, 28, 70 };
		///*
		for (int i = 0; i < value.Length; ++i)
		{
			//invMass[i] = 0;
		}
		//invMass[107] = 0;
		//*/
		for (int i = 0; i < Low_edgeLengths.Length; i++)
		{
			int id0 = Low_edgeIds[2 * i];
			int id1 = Low_edgeIds[2 * i + 1];
			Low_edgeLengths[i] = Vector3.Distance(Low_Xpos[id0], Low_Xpos[id1]);
		}

		for (int i = 0; i < High_edgeLengths.Length; i++)
		{
			int id0 = High_edgeIds[2 * i];
			int id1 = High_edgeIds[2 * i + 1];
			High_edgeLengths[i] = Vector3.Distance(High_Xpos[id0], High_Xpos[id1]);
		}

		graphMesh.vertices = High_Xpos;
		graphMesh.triangles = High_Tris;
		//graphMesh.RecalculateNormals();
		calculateNormal();
		Bounds = graphMesh.bounds;

		meshCollider.sharedMesh = graphMesh;
		//bounds1 = meshCollider.bounds;

	}

	public Vector3[] calculateLow_Normal() {
		Vector3[] normals = new Vector3[Low_numVerts];
		for (int i = 0; i < normals.Length; ++i)
		{
			normals[i] = new Vector3(0.0f, 0.0f, 0.0f);
		}
		for (int i = 0; i < Low_face.Length; i = i + 3)
		{
			Vector3 v1 = Low_Xpos[Low_face[i]];
			Vector3 v2 = Low_Xpos[Low_face[i + 2]];
			Vector3 v3 = Low_Xpos[Low_face[i + 1]];

			Vector3 normal = Vector3.Cross(v2 - v1, v3 - v1).normalized;
			normals[Low_face[i]] += normal;
			normals[Low_face[i + 1]] += normal;
			normals[Low_face[i + 2]] += normal;
		}

		for (int i = 0; i < normals.Length; i++)
		{
			normals[i] = normals[i].normalized;
		}
		return normals;
	}
	public void calculateNormal()
	{
		Vector3[] normals = new Vector3[High_numVerts];
		for (int i = 0; i < normals.Length; ++i)
		{
			normals[i] = new Vector3(0.0f, 0.0f, 0.0f);
		}
		for (int i = 0; i < mesh.tetSurfaceTriIds.Length; i = i + 3)
		{
			Vector3 v1 = High_Xpos[mesh.tetSurfaceTriIds[i]];
			Vector3 v2 = High_Xpos[mesh.tetSurfaceTriIds[i + 2]];
			Vector3 v3 = High_Xpos[mesh.tetSurfaceTriIds[i + 1]];

			Vector3 normal = Vector3.Cross(v2 - v1, v3 - v1).normalized;
			normals[mesh.tetSurfaceTriIds[i]] += normal;
			normals[mesh.tetSurfaceTriIds[i + 1]] += normal;
			normals[mesh.tetSurfaceTriIds[i + 2]] += normal;
		}

		for (int i = 0; i < normals.Length; i++)
		{
			normals[i] = normals[i].normalized;
		}

		graphMesh.normals = normals;
	}

	public float getTetVolume(int tet)
	{
		Matrix3x3 Dm = new Matrix3x3();

		Vector3 X10 = Low_Xpos[Low_Tet[tet * 4 + 1]] - Low_Xpos[Low_Tet[tet * 4 + 0]];
		Vector3 X20 = Low_Xpos[Low_Tet[tet * 4 + 2]] - Low_Xpos[Low_Tet[tet * 4 + 0]];
		Vector3 X30 = Low_Xpos[Low_Tet[tet * 4 + 3]] - Low_Xpos[Low_Tet[tet * 4 + 0]];

		Dm.M11 = X10.x;
		Dm.M12 = X20.x;
		Dm.M13 = X30.x;
		Dm.M21 = X10.y;
		Dm.M22 = X20.y;
		Dm.M23 = X30.y;
		Dm.M31 = X10.z;
		Dm.M32 = X20.z;
		Dm.M33 = X30.z;

		float det = Matrix3x3.Determinant(Dm);
		float volume = det * 0.1666667f; // 体积
		return volume;
	}

	public void preSolve(float dt, Vector3 gravity)
	{
		for (int i = 0; i < Low_numVerts; ++i)
		{
			if (invMass[i] == 0.0)
			{
				continue;
			}
			vel[i] = new Vector3(vel[i].x * damping, vel[i].y * damping, vel[i].z * damping);
			vel[i].x += gravity.x * dt; vel[i].y += gravity.y * dt; vel[i].z += gravity.z * dt;
			Low_prevX[i] = Low_Xpos[i];
			Low_Xpos[i].x += vel[i].x * dt;
			Low_Xpos[i].y += vel[i].y * dt;
			Low_Xpos[i].z += vel[i].z * dt;
			float x = Low_Xpos[i].x;
			float y = Low_Xpos[i].y;
			float z = Low_Xpos[i].z;

			/*
			if (Math.Abs(x) > 5)
			{
				float _x = x / Math.Abs(x);
				Low_Xpos[i] = Low_prevX[i];
				Low_Xpos[i].x = 5 * _x;
			}
			if (Math.Abs(z) > 5)
			{
				float _z = z / Math.Abs(z);
				Low_Xpos[i] = Low_prevX[i];
				Low_Xpos[i].z = 5 * _z;
			}
			*/
			if (y < -5)
			{
				Low_Xpos[i] = Low_prevX[i];
				Low_Xpos[i].y = -5;
			}
			

		}


	}

	public void solve(float dt)
	{
		//if (openEdgeConstrian)
		solveEdges(edgeCompliance, dt);
		if (openVolumeConstrian)
			solveVolumes(volCompliance, dt);
		//solveCollision();
	}
	public void postSolve(float dt)
	{
		for (int i = 0; i < Low_numVerts; i++)
		{
			if (invMass[i] == 0.0)
				continue;

			vel[i].x = -(Low_prevX[i].x - Low_Xpos[i].x) / dt;
			vel[i].y = -(Low_prevX[i].y - Low_Xpos[i].y) / dt;
			vel[i].z = -(Low_prevX[i].z - Low_Xpos[i].z) / dt;
		}
	}

	public void solveEdges(float compliance, float dt)
	{
		float alpha = compliance / dt / dt;

		for (int i = 0; i < Low_edgeLengths.Length; i++)
		{
			int id0 = Low_edgeIds[2 * i];
			int id1 = Low_edgeIds[2 * i + 1];
			float w0 = invMass[id0];
			float w1 = invMass[id1];
			float w = w0 + w1;
			if (w == 0.0)
				continue;

			Vector3 director = Low_Xpos[id0] - Low_Xpos[id1];
			float len = director.magnitude;
			//方向向量
			director = director / len;

			if (openEdgeConstrian)
			{
				if (len >= 0.01f)
				{
					float k = edgeK;
					float restLen = Low_edgeLengths[i];
					float C = len - restLen;
					if (Math.Abs(C) <= 0.000001) C = 0.0f;
					//if (len > 2 * restLen || len < 0.3 * restLen) k = 1.5f * k;

					float s = -k * C / (w + alpha);

					Low_Xpos[id0] += director * s * w0;
					Low_Xpos[id1] += -director * s * w1;
				}
			}
		}
	}
	public void solveHighEdges()
	{
		for (int i = 0; i < High_edgeLengths.Length; i++)
		{
			int id0 = High_edgeIds[2 * i];
			int id1 = High_edgeIds[2 * i + 1];

			Vector3 director = High_Xpos[id0] - High_Xpos[id1];
			float len = director.magnitude;
			//方向向量
			director = director / len;

			if (len >= 0.01f)
			{
					float k = K;
					float restLen = High_edgeLengths[i];
					float C = len - restLen;
					if (Math.Abs(C) <= 0.000001) C = 0.0f;
				//if (len > 2 * restLen || len < 0.3 * restLen) k = 1.5f * k;

				float s = -k * C;

				High_Xpos[id0] += director * s ;
				High_Xpos[id1] += -director * s ;
			}
		}
	}
	public void solveVolumes(float compliance, float dt)
	{
		float alpha = compliance / dt / dt;

		for (int i = 0; i < Low_numTets; i++)
		{
			float w = 0.0f;

			List<Vector3> grads = new List<Vector3>(4);
			for (int j = 0; j < 4; j++)
			{
				int id0 = Low_Tet[4 * i + volIdOrder[j, 0]];
				int id1 = Low_Tet[4 * i + volIdOrder[j, 1]];
				int id2 = Low_Tet[4 * i + volIdOrder[j, 2]];

				Vector3 temp1 = Low_Xpos[id1] - Low_Xpos[id0];
				Vector3 temp2 = Low_Xpos[id2] - Low_Xpos[id0];
				grads.Add(Vector3.Cross(temp1, temp2) / 6.0f);
				float y = grads[j].magnitude;

				w += invMass[Low_Tet[4 * i + j]] * y * y;

			}
			if (w == 0.0)
				continue;

			float vol = getTetVolume(i);
			float restvol = restVol[i];
			float C = vol - restvol;
			if (Math.Abs(C) < 0.000001f) C = 0;
			float s = -volumeK * C / (w + alpha);

			for (int j = 0; j < 4; j++)
			{
				int id = Low_Tet[4 * i + j];
				Vector3 grad;
				grad = grads[j];
				Low_Xpos[id].x += grad.x * s * invMass[id]; Low_Xpos[id].y += grad.y * s * invMass[id]; Low_Xpos[id].z += grad.z * s * invMass[id];
			}
		}
	}


	private void Update()
	{
		float f = Time.unscaledDeltaTime;
		frame += 1;
		duration += f;
		if (duration >= 2)
		{
			fps = frame / duration;
			frame = 0;
			duration = 0;
		}

		if (Input.GetKeyDown(KeyCode.Space))
		{
			test++;
		}

		StartCoroutine(OnMouseDown());

		
		//低模的软体形变
		float dt = 1.0f / 600.0f;
		for (int step = 0; step < 10; step++)
		{
			preSolve(dt, new Vector3(0, -gravity, 0));
			//UpdateCollision();
			if (pickHigh_Vertex != -1)
			{
				Low_Xpos[pickLow_Vertex] = pickLow_VertexPos;
				//StartCoroutine(OnMouseDown());
			}
			Low_CurNormals = calculateLow_Normal();
			solve(dt);
			postSolve(dt);

			//低模映射到高模

			for (int i = 0; i < Low_numVerts; ++i)
			{
				Dictionary<int, Vector3> innerMap = Map[i];
				Dictionary<int, Vector3> innerSecondMap = Map_Second[i];

				//------测试旋转
				/*
				Dictionary<int, Vector3> innerForRotate_Line = ForRotate_Line[i];
				Quaternion initialRotation;// = Quaternion.LookRotation(ForRotate_Line[i], Low_OriginNormals[i]);
				Quaternion CurRotation;// = Quaternion.LookRotation(, finalN1);
				Matrix4x4 rotationMatrix = new Matrix4x4();
				foreach (KeyValuePair<int, Vector3> x in innerForRotate_Line)
				{
					int other_Vert = x.Key;
					Vector3 vector = x.Value;

					initialRotation = Quaternion.LookRotation(vector, Low_OriginNormals[i]);
					CurRotation = Quaternion.LookRotation(Low_Xpos[i] - Low_Xpos[other_Vert], Low_CurNormals[i]);
					// 计算从初始状态到最终状态的旋转变换矩阵
					Quaternion rotationDifference = Quaternion.Inverse(initialRotation) * CurRotation;
					rotationMatrix = Matrix4x4.Rotate(rotationDifference);
				}
				*/

				Quaternion rotation = Quaternion.FromToRotation(Low_OriginNormals[i], Low_CurNormals[i]);
				foreach (KeyValuePair<int, Vector3> item in innerMap)
				{
					int high_Vertex = item.Key; // 高模顶点的索引
					Vector3 distance = item.Value; // 从低模顶点到高模顶点的距离
												   // 计算低模法线旋转的四元数

					distance = rotation * distance;



					//if (i ==0 ||i == 19)
					//High_Xpos[high_Vertex] = rotationMatrix.MultiplyPoint3x4(High_Xpos[high_Vertex]);



					//if (distance >= 1) Debug.Log(distance);
					//High_Xpos[high_Vertex] += vector * (1 - distance) * (1 - distance) * (1 - distance);
					//High_Xpos[high_Vertex] = Low_Xpos[i] - distance;
					High_Xpos[high_Vertex] = Low_Xpos[i] - (distance + Low_CurNormals[i]);
					if (i == 0 || i == 19 ) Debug.DrawLine(Low_Xpos[i], Low_Xpos[i] - (distance + Low_CurNormals[i]), Color.red);
					if (i == 0 || i == 19) Debug.DrawLine(Low_Xpos[i], Low_Xpos[i] + Low_CurNormals[i], Color.red);
					
				}
				foreach (KeyValuePair<int, Vector3> item in innerSecondMap) {
					int high_Vertex = item.Key; // 高模顶点的索引
					Vector3 distance = item.Value; // 从低模顶点到高模顶点的距离
					distance = rotation * distance;
					//if (i == 0) Debug.DrawLine(Low_Xpos[i], Low_Xpos[i] - (distance + Low_CurNormals[i]), Color.cyan);
				}

				

			}
		}
		//solveHighEdges();


		for (int i = 0; i < Low_edgeIds.Length; i += 2) {
			if (Low_edgeIds[i] == 0 || Low_edgeIds[i + 1] == 0)
				Debug.DrawLine(Low_Xpos[Low_edgeIds[i]], Low_Xpos[Low_edgeIds[i + 1]], Color.yellow);
		}

		if (Input.GetMouseButtonDown(0))
		{
			PickVertex();
		}
		

		graphMesh.vertices = High_Xpos;
		graphMesh.triangles = High_Tris;
		//graphMesh.RecalculateNormals();
		calculateNormal();
		meshCollider.sharedMesh = graphMesh;
		//DrawBoundBoxLine(GetBounds(graphMesh.vertices));
	}

	Bounds GetBounds(Vector3[] vertices)
	{
		Vector3 vMin = new Vector3(float.MaxValue, float.MaxValue, float.MaxValue);
		Vector3 vMax = new Vector3(float.MinValue, float.MinValue, float.MinValue);

		for (int i = 0; i < vertices.Length; ++i)
		{
			Vector3 P = vertices[i];
			vMin = Vector3.Min(vMin, P);
			vMax = Vector3.Max(vMax, P);
		}

		Bounds.center = 0.5f * (vMin + vMax);
		Bounds.extents = 0.5f * (vMax - vMin);
		bounds1.center = 0.5f * (vMin + vMax);
		bounds1.extents = 0.5f * (vMax - vMin);
		return bounds1;
	}
	public static void DrawBoundBoxLine(Bounds bounds, Color color = default(Color), float offsetSize = 1f, float duration = 0.1f)
	{
		color = Color.red;
		//先计算出包围盒8个点
		Vector3[] points = new Vector3[8];
		var width_x = bounds.size.x * offsetSize;
		var hight_y = bounds.size.y * offsetSize;
		var length_z = bounds.size.z * offsetSize;

		var LeftBottomPoint = bounds.min;
		var rightUpPoint = bounds.max;
		var centerPoint = bounds.center;
		var topPoint = new Vector3(centerPoint.x, centerPoint.y + hight_y / 2, centerPoint.z);
		var bottomPoint = new Vector3(centerPoint.x, centerPoint.y - hight_y * 0.5f, centerPoint.z);

		points[0] = LeftBottomPoint + Vector3.right * width_x;
		points[1] = LeftBottomPoint + Vector3.up * hight_y;
		points[2] = LeftBottomPoint + Vector3.forward * length_z;

		points[3] = rightUpPoint - Vector3.right * width_x;
		points[4] = rightUpPoint - Vector3.up * hight_y;
		points[5] = rightUpPoint - Vector3.forward * length_z;

		points[6] = LeftBottomPoint;
		points[7] = rightUpPoint;

		Debug.DrawLine(LeftBottomPoint, points[0], color, duration);
		Debug.DrawLine(LeftBottomPoint, points[1], color, duration);
		Debug.DrawLine(LeftBottomPoint, points[2], color, duration);

		Debug.DrawLine(rightUpPoint, points[3], color, duration);
		Debug.DrawLine(rightUpPoint, points[4], color, duration);
		Debug.DrawLine(rightUpPoint, points[5], color, duration);

		Debug.DrawLine(points[1], points[3], color, duration);
		Debug.DrawLine(points[2], points[4], color, duration);
		Debug.DrawLine(points[0], points[5], color, duration);

		Debug.DrawLine(points[2], points[3], color, duration);
		Debug.DrawLine(points[0], points[4], color, duration);
		Debug.DrawLine(points[1], points[5], color, duration);
		//Debug.Log(bounds.extents);

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
			if (collider.gameObject.CompareTag("organ"))
			{
				//Debug.Log(hit.collider.gameObject.name);
				UnityEngine.Mesh mesh0 = collider.sharedMesh;
				Vector3[] vertices = mesh0.vertices;
				int[] triangles = mesh0.triangles;
				Vector3 p0 = hit.transform.TransformPoint(vertices[triangles[hit.triangleIndex * 3]]);
				Vector3 p1 = hit.transform.TransformPoint(vertices[triangles[hit.triangleIndex * 3 + 1]]);
				Vector3 p2 = hit.transform.TransformPoint(vertices[triangles[hit.triangleIndex * 3 + 2]]);
				pickHigh_Vertex = MinDistanceMouse(p0, p1, p2, triangles, hit.triangleIndex);


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


	IEnumerator OnMouseDown()
	{
		if (pickHigh_Vertex != -1)
		{
			Vector3 targetScreenPos = Camera.main.WorldToScreenPoint(High_Xpos[pickHigh_Vertex]);

			
			Vector3 director;
			pickLow_Vertex = mapForPick[pickHigh_Vertex]; //对应的低模的顶点
			Dictionary<int, Vector3> innerDict = Map[pickLow_Vertex]; 
			director = innerDict[pickHigh_Vertex];
			//director = Low_Xpos[pickLow_Vertex] - High_Xpos[pickHigh_Vertex];
			Quaternion rotation = Quaternion.FromToRotation(Low_OriginNormals[pickLow_Vertex], Low_CurNormals[pickLow_Vertex]);
			director = rotation * director;

			while (Input.GetMouseButton(0))
			{
				Vector3 mousePos = new Vector3(Input.mousePosition.x, Input.mousePosition.y, targetScreenPos.z);
				var targetPos = Camera.main.ScreenToWorldPoint(mousePos);
				pickHigh_VertexPos = targetPos;
				pickLow_VertexPos = pickHigh_VertexPos + (director+Low_CurNormals[pickLow_Vertex]);
				yield return new WaitForFixedUpdate();//循环执行
			}
		}
		pickHigh_Vertex = -1;
		pickLow_Vertex = -1;
	}


	private void OnDrawGizmos()
	{
		for (int i = 0; i < Low_numVerts; ++i)
		{
			Gizmos.color = Color.black;
			Gizmos.DrawSphere(Low_Xpos[i], 0.07f);
		}

		if (pickHigh_Vertex != -1)
		{
			Gizmos.color = Color.red;
			Gizmos.DrawSphere(High_Xpos[pickHigh_Vertex], 0.05f);
		}
		
	}
}

