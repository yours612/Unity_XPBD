using UnityEngine;
using System.Collections.Generic;
namespace YanCheZuo
{
    [RequireComponent(typeof(MeshFilter))]
    public class TestMeshViewer : MonoBehaviour
    {
        private MeshFilter _meshFilter;
        private UnityEngine.Mesh _mesh;

        public List<Vector3> verticesList = new List<Vector3>();
        //public List<Vector2> uvList = new List<Vector2>();
        public List<int> triList = new List<int>();

        private void Start()
        {
            _meshFilter = this.GetComponent<MeshFilter>();
            _mesh = _meshFilter.mesh;

            ReadMeshInfo();
        }
        private void ReadMeshInfo()
        {

            for (int i = 0, imax = _mesh.vertexCount; i < imax; ++i)
            {
                verticesList.Add(_mesh.vertices[i]);
                //uvList.Add(_mesh.uv[i]);
            }

            for (int i = 0, imax = _mesh.triangles.Length; i < imax; ++i)
            {
                triList.Add(_mesh.triangles[i]);
            }
        }
    }
}
