using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public struct _Matrix3x3_
{
	public float M11, M12, M13;
	public float M21, M22, M23;
	public float M31, M32, M33;

	public static float Determinant(_Matrix3x3_ matrix)
	{
		return matrix.M11 * matrix.M22 * matrix.M33 +
			   matrix.M12 * matrix.M23 * matrix.M31 +
			   matrix.M13 * matrix.M21 * matrix.M32 -
			   matrix.M31 * matrix.M22 * matrix.M13 -
			   matrix.M32 * matrix.M23 * matrix.M11 -
			   matrix.M33 * matrix.M21 * matrix.M12;
	}
}
public class Low2HighXPBD : MonoBehaviour
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


	//-----------与肝脏接触面--------------------
	OrganContact organContact;

	//-----------低模：四面体tet、顶点vert、边表edge      高模：顶点vert， 三角面tri-----------------

	//低模对应高模的顶点的连线向量----高模的顶点与最近的低模点相连
	Dictionary<int, Dictionary<int, Vector3>> Map = new Dictionary<int, Dictionary<int, Vector3>>();
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

	Vector3[] low_faceOriginNormals;
	Vector3[] low_faceCurNormals;
	Dictionary<int, Dictionary<int, float[]>> FaceToVert = new Dictionary<int, Dictionary<int, float[]>>();
	// （三角面序号 、（高模顶点、投影点数据WA,WB,WC,点到面的距离-有正负））

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
	public Vector3[] High_Xpos;
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


	public void initialize() {
		mesh = new Mesh();
		mesh.GetInfor();

		High_numVerts = mesh.vertsPos.Length / 3;
		High_Xpos = new Vector3[High_numVerts];
		High_X = mesh.vertsPos;
		for (int i = 0, j = 0; i < High_X.Length; i += 3)
		{
			Vector3 vertex;
			vertex.x = High_X[i];
			vertex.y = High_X[i + 1];
			vertex.z = High_X[i + 2];
			High_Xpos[j++] = vertex;
		}
	}
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
															   
		low_faceOriginNormals = new Vector3[Low_face.Length / 3];
		low_faceCurNormals = new Vector3[Low_face.Length / 3];

		//初始化低模高模的顶点数组；
		for (int i = 0, j = 0; i < Low_X.Length; i += 3)
		{
			Vector3 vertex;
			vertex.x = Low_X[i];
			vertex.y = Low_X[i + 1];
			vertex.z = Low_X[i + 2];
			Low_Xpos[j++] = vertex;
			//VirtualXpos[j++] = vertex;
		}
		Low_prevX = new Vector3[Low_Xpos.Length];
		for (int i = 0, j = 0; i < High_X.Length; i += 3)
		{
			Vector3 vertex;
			vertex.x = High_X[i];
			vertex.y = High_X[i + 1];
			vertex.z = High_X[i + 2];
			High_Xpos[j++] = vertex;
		}

		//初始化低模顶点法线\面法线---------------------------------
		Low_OriginNormals = calculateLow_Normal();
		low_faceOriginNormals = calculateLow_FaceNormals();
		//初始化低模对应高模的顶点连线
		for (int i = 0; i < Low_numVerts; ++i)
		{
			Map[i] = new Dictionary<int, Vector3>();
		}
		for (int j = 0; j < High_numVerts; ++j)
		{
			Vector3 pos = High_Xpos[j];
			float minDistance = float.MaxValue;
			int closestLowPolyVertexIndex = -1;

			for (int i = 0; i < Low_numVerts; ++i)
			{
				float distance = Vector3.Distance(pos, Low_Xpos[i]);
				if (distance < minDistance)
				{
					minDistance = distance;
					closestLowPolyVertexIndex = i;
				}
			}

			if (closestLowPolyVertexIndex != -1)
			{
				mapForPick[j] = closestLowPolyVertexIndex;
				Map[closestLowPolyVertexIndex].Add(j, Low_Xpos[closestLowPolyVertexIndex] - pos - Low_OriginNormals[closestLowPolyVertexIndex]);
			}
		}

		//初始化低模三角面控制哪些高模的顶点---------------------------------------
		for (int i = 0; i < Low_face.Length / 3; ++i) { 
			FaceToVert[i] = new Dictionary<int, float[]>();
		}
		for (int j = 0; j < High_numVerts; ++j)
		{
			float minDistance = float.MaxValue;
			Vector3 closestVector = Vector3.zero;
			int closestTriangleIndex = -1;
			float[] W = new float[4];

			Vector3 v1 = new Vector3(), v2 = new Vector3(), v3 = new Vector3(), normal = new Vector3();
			float distance = 0;

			for (int i = 0, index = 0; i < Low_face.Length; i += 3) {
				Vector3 _v1 = Low_Xpos[Low_face[i]];
				Vector3 _v2 = Low_Xpos[Low_face[i + 2]];
				Vector3 _v3 = Low_Xpos[Low_face[i + 1]];
				normal = Vector3.Cross(_v2 - _v1, _v3 - _v1).normalized;

				float _distance = Vector3.Dot(normal, High_Xpos[j] - _v1);
				Vector3 p = High_Xpos[j] - _distance * normal;

				if (IsPointInTriangle(p, _v1, _v2, _v3))//投影在三角面内
				{
					if (Mathf.Abs(_distance) < minDistance)
					{
						minDistance = Mathf.Abs(_distance);
						v1 = _v1; v2 = _v2; v3 = _v3;
						distance = _distance;
						closestTriangleIndex = index;
					}
				}
				else {
					Vector3 closestPoint = ClosestPointToTriangleVerticesAndEdges(High_Xpos[j], _v1, _v2, _v3);
					float outTri_distance = (High_Xpos[j] - closestPoint).magnitude;
					if (Mathf.Abs(outTri_distance) < minDistance)
					{
						minDistance = Mathf.Abs(outTri_distance);
						v1 = _v1; v2 = _v2; v3 = _v3;
						distance = _distance;
						closestTriangleIndex = index;
					}
				}
				index++;
			}
			normal = Vector3.Cross(v2 - v1, v3 - v1).normalized;
			Vector3 toPlane = distance * normal;
			closestVector = toPlane;
			Vector3 projectedPoint = High_Xpos[j] - distance * normal;

			Vector3 X0 = v2 - v1, X1 = v3 - v1, X2 = projectedPoint - v1;
			float d00 = Vector3.Dot(X0, X0);
			float d01 = Vector3.Dot(X0, X1);
			float d11 = Vector3.Dot(X1, X1);
			float d20 = Vector3.Dot(X2, X0);
			float d21 = Vector3.Dot(X2, X1);
			float denom = d00 * d11 - d01 * d01;
			float wB = (d11 * d20 - d01 * d21) / denom;
			float wC = (d00 * d21 - d01 * d20) / denom;
			float wA = 1.0f - wC - wB;
			W = new float[] { wA, wB, wC, distance };

			FaceToVert[closestTriangleIndex].Add(j, W);
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

		int[] value = { 5, 3, 7, 9, 15 ,23 ,27};
		///*
		for (int i = 0; i < value.Length; ++i)
		{
			//invMass[i] = 0;
		}
		//invMass[25] = 0;

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


		//--------接触面
		organContact = new OrganContact();
		organContact.initialize();
		//organContact
	}

	public Vector3[] calculateLow_FaceNormals()
	{
		Vector3[] normals = new Vector3[Low_face.Length/3];
		for (int i = 0; i < normals.Length; ++i)
		{
			normals[i] = new Vector3(0.0f, 0.0f, 0.0f);
		}
		for (int i = 0; i < Low_face.Length / 3; i++)
		{
			Vector3 v1 = Low_Xpos[Low_face[i*3]];
			Vector3 v2 = Low_Xpos[Low_face[i*3 + 2]];
			Vector3 v3 = Low_Xpos[Low_face[i*3 + 1]];

			Vector3 normal = Vector3.Cross(v2 - v1, v3 - v1).normalized;
			normals[i] = normal.normalized;
		}
		return normals;
	}

	public Vector3[] calculateLow_Normal()
	{
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
		
			if (y < -5)
			{
				Low_Xpos[i] = Low_prevX[i];
				Low_Xpos[i].y = -5;
			}
			*/

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


		//接触面先更新
		for (int i = 0; i < organContact.num; ++i)
		{
			High_Xpos[organContact.contact_gallbladder[i]] = organContact.contact_Pos[i];
		}

		//低模的软体形变
		float dt = 1.0f / 600.0f;
		for (int step = 0; step < 10; step++)
		{
			preSolve(dt, new Vector3(0, -gravity, 0));
			//UpdateCollision();
			if (pickHigh_Vertex != -1)
			{
				Low_Xpos[pickLow_Vertex] = pickLow_VertexPos;
			}
			Low_CurNormals = calculateLow_Normal();
			solve(dt);
			postSolve(dt);

			//低模映射到高模
			for (int i = 0; i < Low_face.Length/3; i++) {
				Dictionary<int, float[]> inner = FaceToVert[i];
				Vector3 _v1 = Low_Xpos[Low_face[i*3]];
				Vector3 _v2 = Low_Xpos[Low_face[i*3 + 2]];
				Vector3 _v3 = Low_Xpos[Low_face[i*3 + 1]];
				Vector3 normal = Vector3.Cross(_v2 - _v1, _v3 - _v1).normalized;
				foreach (KeyValuePair<int, float[]> item in inner)
				{
					int high_Vertex = item.Key; // 高模顶点的索引
					float[] W = item.Value; // 从低模顶点到高模顶点的距离
	
					Vector3 p = _v1 * W[0] + _v2 * W[1] + _v3 * W[2];

					High_Xpos[high_Vertex] = p + normal * W[3];
				}
			}

			//与肝脏的接触面的点
			for (int i = 0; i < organContact.num; ++i)
			{
				organContact.contact_Pos[i] = High_Xpos[organContact.contact_gallbladder[i]];
			}

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
				//pickLow_VertexPos = pickHigh_VertexPos + (director + Low_CurNormals[pickLow_Vertex]);
				pickLow_VertexPos = Low_Xpos[pickLow_Vertex] + (targetPos - High_Xpos[pickHigh_Vertex]);
				yield return new WaitForFixedUpdate();//循环执行
			}
		}
		pickHigh_Vertex = -1;
		pickLow_Vertex = -1;
	}


	private void OnDrawGizmos()
	{
		/*
		int[] point = { 392, 123, 124,236,309,25,126 };
		for (int i = 0; i < point.Length; ++i)
		{

				Gizmos.color = Color.green;
				Gizmos.DrawSphere(High_Xpos[point[i]], 0.15f);
			
		}
		*/

		if (pickHigh_Vertex != -1)
		{
			Gizmos.color = Color.red;
			Gizmos.DrawSphere(High_Xpos[pickHigh_Vertex], 0.05f);
		}

	}

	Vector3 ClosestPointToTriangleVerticesAndEdges(Vector3 P, Vector3 A, Vector3 B, Vector3 C)
	{
		// 初始化最近点为A，这只是一个起始比较值
		Vector3 closestPoint = A;
		float minDistanceSquared = (P - A).sqrMagnitude;

		// 检查B和C是否更近
		float distanceSquaredB = (P - B).sqrMagnitude;
		if (distanceSquaredB < minDistanceSquared)
		{
			closestPoint = B;
			minDistanceSquared = distanceSquaredB;
		}

		float distanceSquaredC = (P - C).sqrMagnitude;
		if (distanceSquaredC < minDistanceSquared)
		{
			closestPoint = C;
			minDistanceSquared = distanceSquaredC;
		}

		// 检查点到每条边的最近点
		Vector3[] edgeStartPoints = { A, B, C };
		Vector3[] edgeEndPoints = { B, C, A };
		for (int i = 0; i < 3; i++)
		{
			Vector3 closestPointOnEdge = ClosestPointOnLineSegment(edgeStartPoints[i], edgeEndPoints[i], P);
			float distanceSquaredEdge = (P - closestPointOnEdge).sqrMagnitude;
			if (distanceSquaredEdge < minDistanceSquared)
			{
				closestPoint = closestPointOnEdge;
				minDistanceSquared = distanceSquaredEdge;
			}
		}

		return closestPoint;
	}

	bool IsPointInTriangle(Vector3 P, Vector3 A, Vector3 B, Vector3 C)
	{
		Vector3 v0 = C - A;
		Vector3 v1 = B - A;
		Vector3 v2 = P - A;

		float d00 = Vector3.Dot(v0, v0);
		float d01 = Vector3.Dot(v0, v1);
		float d11 = Vector3.Dot(v1, v1);
		float d20 = Vector3.Dot(v2, v0);
		float d21 = Vector3.Dot(v2, v1);
		float denom = d00 * d11 - d01 * d01;

		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0f - v - w;

		return u >= 0 && v >= 0 && w >= 0;
	}
	Vector3 ClosestPointOnLineSegment(Vector3 A, Vector3 B, Vector3 P)
	{
		Vector3 AP = P - A;
		Vector3 AB = B - A;
		float magnitudeAB = AB.sqrMagnitude;  // 长度的平方
		float ABAPproduct = Vector3.Dot(AP, AB);
		float distance = ABAPproduct / magnitudeAB;

		if (distance < 0)
		{
			return A;
		}
		else if (distance > 1)
		{
			return B;
		}
		else
		{
			return A + AB * distance;
		}
	}

}