using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public struct Matrix3x3
{
	public float M11, M12, M13;
	public float M21, M22, M23;
	public float M31, M32, M33;

	public static float Determinant(Matrix3x3 matrix)
	{
		return matrix.M11 * matrix.M22 * matrix.M33 +
			   matrix.M12 * matrix.M23 * matrix.M31 +
			   matrix.M13 * matrix.M21 * matrix.M32 -
			   matrix.M31 * matrix.M22 * matrix.M13 -
			   matrix.M32 * matrix.M23 * matrix.M11 -
			   matrix.M33 * matrix.M21 * matrix.M12;
	}
}

public class XPBD : MonoBehaviour
{
	int x = 1;

	public float K = 0.001f;
	public bool openEdgeConstrian = true;
	public float edgeK = 0.001f;
	public bool openVolumeConstrian = false;
	public float volumeK = 1.0f;
	public bool openVirtualPosdEdgeConstrian = true;
	public float virtualK = 0.05f;
	public bool openEdgeAngleConstrian = true;
	public float edgeAngleK = 0.05f;
	public bool openOppositeEdgeConstrian = true;
	public float oppositeEdgeK = 0.05f;
	public float edgeCompliance;
	public float volCompliance;
	public float gravity;
	public float damping = 0.99f;

	//float virtualK = 0.1f;

	int numVerts;
	int numTets;
	float[] X;
	Vector3[] Xpos;
	Vector3[] prevX;
	Vector3[] vel;
	Vector3[] VirtualXpos;
	int[] Tris;

	int[] Tet;
	int[] edgeIds;
	float[] restVol;
	float[] restSurVol;
	float[] edgeLengths;
	float[] oppositeedgeLengths;
	float[] virtualEdgeLengths;
	float[] invMass;	

	public Mesh mesh;
	//bunnyMesh Mesh;
	float meshBallR = 0.1f;

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

	public int pickVertex = -1;
	Vector3 pickVertexPos = new Vector3(0, 0, 0);

	//fps
	public float fps;
	int frame;
	float duration;

	//法线约束-测试
	Vector3[] Normal;
	//表面体积约束（整体体积）
	float restSurfaceVolume;//初始体积
	float stiffness;
	int[][] map; //map是每个顶点所在三角面的序列。如map[0]是点0所在的是第几个三角面。
	Vector3[] taop; //顶点法线


	// 初始化参数
	void Start()
	{
		graphMesh = GetComponent<MeshFilter>().mesh;
		meshCollider = GetComponent<MeshCollider>();

		mesh = new Mesh();
		mesh.GetInfor();

		numVerts = mesh.vertsPos.Length / 3;
		numTets = mesh.tetId.Length / 4;

		X = mesh.vertsPos;
		Xpos = new Vector3[numVerts];
		VirtualXpos = new Vector3[numVerts];
		vel = new Vector3[numVerts];
		Tet = mesh.tetId;
		edgeIds = mesh.tetEdgeIds;
		restVol = new float[numTets];
		invMass = new float[numVerts];
		//restSurVol = new float[mesh.tetSurfaceTriIds.Length / 3];
		edgeLengths = new float[edgeIds.Length / 2];
		//oppositeedgeLengths = new float[mesh.oppositeEdgesId.Length / 2];
		Normal = new Vector3[numVerts];
		virtualEdgeLengths = new float[edgeIds.Length];

		//edgeCompliance = 0;
		//volCompliance = 0;

		for (int i = 0, j = 0; i < X.Length; i += 3)
		{
			Vector3 vertex;
			vertex.x = X[i] * 0.5f;
			vertex.y = X[i + 1] * 0.5f;
			vertex.z = X[i + 2] * 0.5f;
			Xpos[j] = vertex;
			VirtualXpos[j++] = vertex;
		}
		prevX = new Vector3[Xpos.Length];

		Tris = new int[mesh.tetSurfaceTriIds.Length];
		for (int i = 0; i < mesh.tetSurfaceTriIds.Length; i += 3) {
			Tris[i] = mesh.tetSurfaceTriIds[i];
			Tris[i+1] = mesh.tetSurfaceTriIds[i+2];
			Tris[i+2] = mesh.tetSurfaceTriIds[i+1];
		}

		//所有数组置零
		for (int i = 0; i < numVerts; ++i) {
			//VirtualXpos[i] = new Vector3(0.0f, 0.0f, 0.0f);
			vel[i] = new Vector3(0.0f, 0.0f, 0.0f);
			restVol[i] = 0;
			invMass[i] = 0;
		}

		//初始化质量矩阵，剩余体积矩阵, 边长矩阵
		for (int i = 0; i < numTets; i++)
		{
			float vol = getTetVolume(i);
			restVol[i] = vol;
			float pInvMass = vol > 0.0f ? 1.0f / (vol / 4.0f) : 0.0f;
			invMass[Tet[4 * i]] += pInvMass;
			invMass[Tet[4 * i + 1]] += pInvMass;
			invMass[Tet[4 * i + 2]] += pInvMass;
			invMass[Tet[4 * i + 3]] += pInvMass;
		}

		int[] value = { 1, 0, 8, 84, 28, 70 };
		///*
		for (int i = 0; i < value.Length; ++i)
		{
			//invMass[i] = 0;
		}
		invMass[15] = 0;
		//*/
		for (int i = 0; i < edgeLengths.Length; i++)
		{
			int id0 = edgeIds[2 * i];
			int id1 = edgeIds[2 * i + 1];
			edgeLengths[i] = Vector3.Distance(Xpos[id0], Xpos[id1]);
		}

		/*
		for (int i = 0; i < oppositeedgeLengths.Length; ++i) {
			int id0 = mesh.oppositeEdgesId[2 * i];
			int id1 = mesh.oppositeEdgesId[2 * i + 1];
			oppositeedgeLengths[i] = Vector3.Distance(Xpos[id0], Xpos[id1]);
		}
		*/

		graphMesh.vertices = Xpos;
		graphMesh.triangles = Tris;
		//graphMesh.RecalculateNormals();
		calculateNormal();
		Bounds = graphMesh.bounds;

		meshCollider.sharedMesh = graphMesh;
		bounds1 = meshCollider.bounds;

		//初始法线方向、虚拟节点构建
		
		for (int i = 0; i < numVerts; ++i) {
			Normal[i] = graphMesh.normals[i];

			float ratio = -0.5f;
			Vector3 v = Xpos[i];
			v += Normal[i] * ratio;

			VirtualXpos[i] = v;
		}

		for (int i = 0; i < edgeLengths.Length; i++) {
			int id0 = edgeIds[2 * i];
			int id1 = edgeIds[2 * i + 1];
			virtualEdgeLengths[i*2] = Vector3.Distance(Xpos[id0], VirtualXpos[id1]);
			virtualEdgeLengths[i * 2+1] = Vector3.Distance(Xpos[id1], VirtualXpos[id0]);
		}


			//表面体积约束
			restSurfaceVolume = 0;		
		for (int i = 0; i < Tris.Length / 3; ++i) {
			Vector3 p1 = Xpos[Tris[3 * i]];
			Vector3 p2 = Xpos[Tris[3 * i+1]];
			Vector3 p3 = Xpos[Tris[3 * i+2]];
			restSurfaceVolume += triangleVolume(p1, p2, p3);
		}
		map = new int[numVerts][];
		for (int i = 0; i < numVerts; ++i) {
			List<int> index = new List<int>();

			for (int j = 0; j < Tris.Length / 3; j++) {
				if (Tris[3 * j] == i || Tris[3 * j + 1] == i || Tris[3 * j + 2] == i)
				{
					index.Add(j);
				}
			}
			map[i] = index.ToArray();
		}

		
	}

	public void calculateNormal() {
		int j = 0;
		Vector3[] normals = new Vector3[numVerts];
		for (int i = 0; i < normals.Length; ++i) {
			normals[i] = new Vector3(0.0f, 0.0f, 0.0f);
		}
		for (int i = 0; i < mesh.tetSurfaceTriIds.Length; i = i+3) {
			Vector3 v1 = Xpos[mesh.tetSurfaceTriIds[i]];
			Vector3 v2 = Xpos[mesh.tetSurfaceTriIds[i+2]];
			Vector3 v3 = Xpos[mesh.tetSurfaceTriIds[i+1]];

			Vector3 normal = Vector3.Cross(v2 - v1, v3 - v1).normalized;
			normals[mesh.tetSurfaceTriIds[i]] += normal;
			normals[mesh.tetSurfaceTriIds[i+1]] += normal;
			normals[mesh.tetSurfaceTriIds[i+2]] += normal;
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

		Vector3 X10 = Xpos[Tet[tet * 4 + 1]] - Xpos[Tet[tet * 4 + 0]];
		Vector3 X20 = Xpos[Tet[tet * 4 + 2]] - Xpos[Tet[tet * 4 + 0]];
		Vector3 X30 = Xpos[Tet[tet * 4 + 3]] - Xpos[Tet[tet * 4 + 0]];

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

	public void preSolve(float dt, Vector3 gravity) {
		for (int i = 0; i < numVerts; ++i)
		{
			if (invMass[i] == 0.0)
			{
				continue;
			}
			vel[i] = new Vector3( vel[i].x * damping, vel[i].y * damping, vel[i].z * damping);
			vel[i].x += gravity.x * dt; vel[i].y += gravity.y * dt; vel[i].z += gravity.z * dt;
			prevX[i] = Xpos[i];
			Xpos[i].x += vel[i].x * dt;
			Xpos[i].y += vel[i].y * dt;
			Xpos[i].z += vel[i].z * dt;
			float x = Xpos[i].x;
			float y = Xpos[i].y;
			float z = Xpos[i].z;

			if (Math.Abs(x) > 5)
			{
				float _x = x / Math.Abs(x);
				Xpos[i] = prevX[i];
				Xpos[i].x = 5 * _x;
			}
			if (Math.Abs(z) > 5)
			{
				float _z = z / Math.Abs(z);
				Xpos[i] = prevX[i];
				Xpos[i].z = 5 * _z;
			}
			if (y < -2)
			{
				Xpos[i] = prevX[i];
				Xpos[i].y = -2;
			}

		}


	}

	public void solve(float dt)
	{
		//if (openEdgeConstrian)
		solveEdges(edgeCompliance, dt);
		if (openOppositeEdgeConstrian)
			solveOppositeEdges(dt);
		if (openVolumeConstrian)
			solveVolumes(volCompliance, dt);
		//solveNormals(dt);
		//solveCollision();
		//solveSurVolume();
	}
	public void postSolve(float dt)
		{
			for (int i = 0; i < numVerts; i++)
			{
				if (invMass[i] == 0.0)
					continue;

				vel[i].x = -(prevX[i].x - Xpos[i].x) / dt;
				vel[i].y = -(prevX[i].y - Xpos[i].y) / dt;
				vel[i].z = -(prevX[i].z - Xpos[i].z) / dt;
			}
		}

	public void solveEdges(float compliance, float dt)
	{
		float alpha = compliance / dt / dt;
		
		for (int i = 0; i < edgeLengths.Length; i++)
		{
			int id0 = edgeIds[2 * i];
			int id1 = edgeIds[2 * i + 1];
			float w0 = invMass[id0];
			float w1 = invMass[id1];
			float w = w0 + w1;
			if (w == 0.0)
				continue;
			
			Vector3 director = Xpos[id0] - Xpos[id1];
			float len = director.magnitude;
			//方向向量
			director = director / len;

			if (openEdgeConstrian)
			{
				if (len >= 0.01f)
				{
					float k = edgeK;
					float restLen = edgeLengths[i];
					float C = len - restLen;
					if (Math.Abs(C) <= 0.000001) C = 0.0f;
					//if (len > 2 * restLen || len < 0.3 * restLen) k = 1.5f * k;

					float s = - k * C / (w + alpha);

					Xpos[id0] += director * s * w0;
					Xpos[id1] += -director * s * w1;
				}
			}



			Vector3 id0director = Xpos[id0] - VirtualXpos[id0];
			Vector3 id1director = Xpos[id1] - VirtualXpos[id1];
			float len0 = id0director.magnitude;
			float len1 = id1director.magnitude;
			id0director = id0director / len0;
			id1director = id1director / len1;

			Vector3 id01director = Xpos[id0] - VirtualXpos[id1];
			Vector3 id10director = Xpos[id1] - VirtualXpos[id0];
			float len01 = id01director.magnitude;
			float len10 = id10director.magnitude;
			id01director /= len01;
			id10director /= len10;

			if (id0 == 183) {
				Debug.DrawLine(VirtualXpos[id0], Xpos[id0], Color.red);
				Debug.DrawLine(VirtualXpos[id0], VirtualXpos[id0] + Normal[id0], Color.blue);
			}

			if (openVirtualPosdEdgeConstrian)
			{
				//虚拟节点--------------自身与自身
				//Vector3 id0director = Xpos[id0] - VirtualXpos[id0]; //{ Xpos[id0].x - VirtualXpos[id0].x, Xpos[id0].y - VirtualXpos[id0].y, Xpos[id0].z - VirtualXpos[id0].z };
				//Vector3 id1director = Xpos[id1] - VirtualXpos[id1];//{ Xpos[id1].x - VirtualXpos[id1].x, Xpos[id1].y - VirtualXpos[id1].y, Xpos[id1].z - VirtualXpos[id1].z };
				//float len0 = id0director.magnitude;// XMVectorGetX(XMVector3Length(XMLoadFloat3(&Xpos[id0]) - XMLoadFloat3(&VirtualXpos[id0])));
				//float len1 = id1director.magnitude; //XMVectorGetX(XMVector3Length(XMLoadFloat3(&Xpos[id1]) - XMLoadFloat3(&VirtualXpos[id1])));

				//id0director = id0director / len0;
				//id1director = id1director / len1;
				if (Math.Abs(len0) >= 0.00001)
				{
					float restLen = 0.5f;
					float C = len0 - restLen;
					if (Math.Abs(C) >= 0.000001) {
						float k = virtualK;
						float s = -k * C / (1 * w0);
						Xpos[id0] += id0director * s * w0;
					}
				}
				if (Math.Abs(len1) >= 0.00001)
				{
					float restLen = 0.5f;
					float C = len1 - restLen;
					if (Math.Abs(C) >= 0.000001) {
						float k = virtualK;
						float s = -k * C / (1 * w1);
						Xpos[id1] += id1director * s * w1;
					}
				}
				//虚拟节点--------------自身与另外
				if (Math.Abs(len01) >= 0.00001) {
					float restlen = virtualEdgeLengths[i*2];
					float C = len01 - restlen;
					if (Math.Abs(C) >= 0.000001) {
						float k = virtualK;
						float s = -k * C / (1 * w0);
						Xpos[id0] += id01director * s * w0;
					}
				}
				if (Math.Abs(len10) >= 0.00001)
				{
					float restlen = virtualEdgeLengths[i * 2+1];
					float C = len10 - restlen;
					if (Math.Abs(C) >= 0.000001)
					{
						float k = virtualK;
						float s = -k * C / (1 * w1);
						Xpos[id1] += id10director * s * w1;
					}
				}
			}



			if (openEdgeAngleConstrian) { 
				///*
				//边的角度约束
				Vector3 originEdge = VirtualXpos[id0] - VirtualXpos[id1];
				float angle = Vector3.Angle(originEdge, director);
				Vector3 proj;
				if (angle != 0)
				{
					proj = Vector3.Project(originEdge, director);
					Vector3 Tangent = originEdge - proj;

					if (pickVertex != -1)
					{
						if ((id0 == 24 && id1 == 27) || (id0 == 27 && id1 == 24))
						{
							//Debug.DrawLine(VirtualXpos[id0], VirtualXpos[id1], Color.green);
							//Debug.DrawLine(Xpos[id0], Xpos[id1], Color.red);
							//Debug.DrawLine(Xpos[id0], Xpos[id0] - proj, Color.black);
						}
					}

					float C = angle;
					float s = -edgeAngleK * C / (w);
					Xpos[id0] -= Tangent * s * w0;
					Xpos[id1] += Tangent * s * w1;
				}

			}

			/*
			//法线约束
			if (Math.Abs(len0) >= 0.00001)
			{
				Vector3 normal0 = Normal[id0];
				Vector3 Proj0 = Vector3.Project(normal0, id0director);
				Vector3 Tangent0 = normal0 - Proj0;

				float angle0 = Vector3.Angle(normal0, Proj0);
				
				if (id0 == 183)
                {
					Debug.DrawLine(Xpos[id0], Xpos[id0] + Tangent0, Color.yellow);
				}
				
				if (angle0 != 0)
				{
					float k = 0.00001f;
					float C = angle0;
					float s = -k * C / (2 * w0);
					Xpos[id0] -= Tangent0 * s * w0;
				}
			}
			if (Math.Abs(len1) >= 0.00001)
			{
				Vector3 normal1 = Normal[id1];
				Vector3 Proj1 = Vector3.Project(normal1, id1director);
				Vector3 Tangent1 = normal1 - Proj1;

				float angle1 = Vector3.Angle(normal1, Proj1);
				if (angle1 != 0)
				{
					float k = 0.00001f;
					float C = angle1;
					float s = -k * C / (2 * w1);
					Xpos[id1] -= Tangent1 * s * w1;
				}
			}
			*/
			//*/


			/*
			Vector3 id01director = Xpos[id0] - VirtualXpos[id1]; // { Xpos[id0].x - VirtualXpos[id1].x, Xpos[id0].y - VirtualXpos[id1].y, Xpos[id0].z - VirtualXpos[id1].z };
			Vector3 id10director = Xpos[id1] - VirtualXpos[id0];//{ Xpos[id1].x - VirtualXpos[id0].x, Xpos[id1].y - VirtualXpos[id0].y, Xpos[id1].z - VirtualXpos[id0].z };
			float len01 = id01director.magnitude;// XMVectorGetX(XMVector3Length(XMLoadFloat3(&Xpos[id0]) - XMLoadFloat3(&VirtualXpos[id1])));
			float len10 = id10director.magnitude; // XMVectorGetX(XMVector3Length(XMLoadFloat3(&Xpos[id1]) - XMLoadFloat3(&VirtualXpos[id0])));
			id01director = id01director / len01;
			id10director = id10director / len10;
			if (len01 != 0.0)
			{
				float restLen = (float)Math.Sqrt(edgeLengths[i] * edgeLengths[i]);
				float C = len01 - restLen;
				if (Math.Abs(C) <= 0.000001) C = 0.0f;
				float k = virtualK;
				float s = -k * 0.7f * C / (2 * w0);
				Xpos[id0].x += id01director.x * s * w0; Xpos[id0].y += id01director.y * s * w0; Xpos[id0].z += id01director.z * s * w0;
			}
			if (len10 != 0.0)
			{
				float restLen = (float)Math.Sqrt(edgeLengths[i] * edgeLengths[i]);
				float C = len10 - restLen;
				if (Math.Abs(C) <= 0.000001) C = 0.0f;
				float k = virtualK;
				float s = -k * 0.7f * C / (2 * w1);
				Xpos[id1].x += id10director.x * s * w1; Xpos[id1].y += id10director.y * s * w1; Xpos[id1].z += id10director.z * s * w1;
			}
			*/
			
		}
	}
	public void solveVolumes(float compliance, float dt)
	{
		float alpha = compliance / dt / dt;

		for (int i = 0; i < numTets; i++)
		{
			float w = 0.0f;

			List<Vector3> grads = new List<Vector3>(4);
			for (int j = 0; j < 4; j++)
			{
				int id0 = Tet[4 * i + volIdOrder[j,0]];
				int id1 = Tet[4 * i + volIdOrder[j,1]];
				int id2 = Tet[4 * i + volIdOrder[j,2]];

				Vector3 temp1 = Xpos[id1] - Xpos[id0];
				Vector3 temp2 = Xpos[id2] - Xpos[id0];
				grads.Add(Vector3.Cross(temp1, temp2) / 6.0f);
				float y = grads[j].magnitude;

				w += invMass[Tet[4 * i + j]] * y * y;

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
				int id = Tet[4 * i + j];
				Vector3 grad;
				grad = grads[j];
				Xpos[id].x += grad.x * s * invMass[id]; Xpos[id].y += grad.y * s * invMass[id]; Xpos[id].z += grad.z * s * invMass[id];
			}
		}
	}

	public void solveOppositeEdges(float dt) {
		for (int i = 0; i < oppositeedgeLengths.Length; ++i) {
			int id0 = mesh.oppositeEdgesId[2 * i];
			int id1 = mesh.oppositeEdgesId[2 * i + 1];
			float w0 = invMass[id0];
			float w1 = invMass[id1];
			float w = w0 + w1;
			if (w == 0.0)
				continue;

			Vector3 director = Xpos[id0] - Xpos[id1];
			float len = director.magnitude;
			//方向向量
			director = director / len;
			if (len != 0.0)
			{
				float restLen = oppositeedgeLengths[i];
				float C = len - restLen;
				if (Math.Abs(C) <= 0.001) C = 0.0f;
				float s = -edgeK * C / (w);
				Xpos[id0].x += director.x * s * w0; Xpos[id0].y += director.y * s * w0; Xpos[id0].z += director.z * s * w0;
				Xpos[id1].x += -director.x * s * w1; Xpos[id1].y += -director.y * s * w1; Xpos[id1].z += -director.z * s * w1;
			}
		}
	}

	public void solveSurVolume() {
		float cp = 0;
		for (int i = 0; i < Tris.Length / 3; i++)
		{
			Vector3 p1 = Xpos[Tris[3 * i]];
			Vector3 p2 = Xpos[Tris[3 * i + 1]];
			Vector3 p3 = Xpos[Tris[3 * i + 2]];
			cp += triangleVolume(p1, p2, p3); //整体体积
		}
		cp -= K * restSurfaceVolume;

		float sdown = 0;
		for (int i = 0; i < graphMesh.normals.Length; ++i) {
			sdown += invMass[i] * graphMesh.normals[i].sqrMagnitude; 
		}
		float s = cp / sdown;

		for (int i = 0; i < numVerts; ++i) { 
			Xpos[i] -= s * invMass[i] * graphMesh.normals[i];
		}
	}

	void solveNormals(float dt) {
		for (int i = 0; i < numVerts; ++i) {
			Vector3 originN = Normal[i];
			Vector3 nowE = Xpos[i] - VirtualXpos[i];

			if (nowE.magnitude >= 0.00001)
			{
				float angle = Vector3.Angle(originN, nowE);
				Vector3 proj;
				if (angle != 0)
				{
					float w = invMass[i];
					proj = Vector3.Project(originN, nowE);
					Vector3 Tangent = originN - proj;

					if (pickVertex != -1)
					{
						if (i == 381)
						{
							Debug.DrawLine(VirtualXpos[i], VirtualXpos[i] + originN, Color.green);
							Debug.DrawLine(Xpos[i], Xpos[i]+ Tangent, Color.blue);
						}
					}

					float C = angle;
					float s = - K * C / (w);
					Xpos[i] -= Tangent * s * w;
				}
			}
		}


	}
	/*
	 * public float getSurVolume(int i) {
		int id0 = mesh.tetSurfaceTriIds[i * 3];
		int id1 = mesh.tetSurfaceTriIds[i * 3 + 1];
		int id2 = mesh.tetSurfaceTriIds[i * 3 + 2];
		float volume = Vector3.Dot(Vector3.Cross(Xpos[id0], Xpos[id2]), Xpos[id1]);
		return volume;
	}
	 * 
	public void solveBalloons(float dt) {
		Vector3 O = new Vector3(0, 0, 0);

		for (int i = 0; i < mesh.tetSurfaceTriIds.Length / 3; i++) 
		{
			float w = 0.0f;

			List<Vector3> grads = new List<Vector3>(3);
			int id0 = mesh.tetSurfaceTriIds[i*3];
			int id1 = mesh.tetSurfaceTriIds[i*3+1];
			int id2 = mesh.tetSurfaceTriIds[i*3+2];
			//1
			Vector3 temp1 = Xpos[id1] - O;
			Vector3 temp2 = Xpos[id2] - O;
			grads.Add(Vector3.Cross(temp1, temp2) / 6.0f);
			float y = grads[0].magnitude;
			w += invMass[id0] * y * y;
			//2
			temp1 = Xpos[id0] - O;
			temp2 = Xpos[id2] - O;
			grads.Add(Vector3.Cross(temp1, temp2) / 6.0f);
			y = grads[1].magnitude;
			w += invMass[id1] * y * y;
			//3
			temp1 = Xpos[id0] - O;
			temp2 = Xpos[id1] - O;
			grads.Add(Vector3.Cross(temp1, temp2) / 6.0f);
			y = grads[2].magnitude;
			w += invMass[id2] * y * y;

			if (w == 0.0)
				continue;

			float vol = Vector3.Dot( Vector3.Cross(Xpos[id0], Xpos[id2]), Xpos[id1]);
			float restvol = restSurVol[i];
			float C = vol - restvol;
			float s = - C /(w);

			Xpos[id0] += grads[0] * s * invMass[id0];
			Xpos[id1] += grads[1] * s * invMass[id1];
			Xpos[id2] += grads[2] * s * invMass[id2];
		}
	}
	*/

	private float triangleVolume(Vector3 p1, Vector3 p2, Vector3 p3)
	{
		return Vector3.Dot(Vector3.Cross(p1, p2), p3);
	}

	private void Update()
	{
		float f = Time.unscaledDeltaTime;
		frame += 1;
		duration += f;
		if (duration >= 2) {
			fps = frame / duration;
			frame = 0;
			duration = 0;
		}


		StartCoroutine(OnMouseDown());
		/*
		if (Input.GetKeyDown(KeyCode.Space))
		{
			for (int i = 0; i < Xpos.Length; i++)
				Xpos[i].y += 0.5f;
		}
		*/
		//
		
		if (Input.GetKeyDown(KeyCode.Space))
		{
			x++;
		}
		Debug.DrawLine(Xpos[Tet[x*4+0]], Xpos[Tet[x * 4 + 1]], Color.green);
		Debug.DrawLine(Xpos[Tet[x * 4 + 0]], Xpos[Tet[x * 4 + 2]], Color.green);
		Debug.DrawLine(Xpos[Tet[x * 4 + 0]], Xpos[Tet[x * 4 + 3]], Color.green);
		Debug.DrawLine(Xpos[Tet[x * 4 + 3]], Xpos[Tet[x * 4 + 1]], Color.green);
		Debug.DrawLine(Xpos[Tet[x * 4 + 2]], Xpos[Tet[x * 4 + 1]], Color.green);
		Debug.DrawLine(Xpos[Tet[x * 4 + 3]], Xpos[Tet[x * 4 + 2]], Color.green);

		float dt = 1.0f / 600.0f;
		for (int step = 0; step < 10; step++)
		{
			preSolve(dt, new Vector3(0,-gravity,0));
			//UpdateCollision();
			if (pickVertex != -1)
			{
				Xpos[pickVertex] = pickVertexPos;
				//StartCoroutine(OnMouseDown());
			}
			solve(dt);
			postSolve(dt);
		}

		if (Input.GetMouseButtonDown(0)) {
			PickVertex();
		}

		graphMesh.vertices = Xpos;
		graphMesh.triangles = Tris;
		//graphMesh.RecalculateNormals();
		calculateNormal();
		meshCollider.sharedMesh = graphMesh;
		//DrawBoundBoxLine(GetBounds(graphMesh.vertices));
	}

	Bounds GetBounds(Vector3[] vertices) {
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

	void PickVertex() {
		RaycastHit hit;
		Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

		if (Physics.Raycast(ray, out hit))
		{
			if (hit.collider == null)
				return;
			MeshCollider collider = (MeshCollider)hit.collider;
			if (collider.gameObject.CompareTag("organ")) {
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
				//Xpos[pickVertex] = targetPos;
				yield return new WaitForFixedUpdate();//循环执行
			}
		}
		pickVertex = -1;
	}


    private void OnDrawGizmos()
    {
		if (pickVertex != -1)
		{
			Gizmos.color = Color.red;
			Gizmos.DrawSphere(Xpos[pickVertex], 0.05f);
		}
	}
}

