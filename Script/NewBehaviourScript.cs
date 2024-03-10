using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NewBehaviourScript : MonoBehaviour
{
    private Mesh m;
    UnityEngine.Mesh mesh;
    void Start()
    {
        m = new Mesh();
        //Debug.Log(m.GetVPos().Length);
        mesh = GetComponent<MeshFilter>().mesh;
        mesh.vertices = m.GetVPos();
        mesh.triangles = m.tetSurfaceTriIds;
        mesh.RecalculateNormals();
    }

    // Update is called once per frame
    void Update()
    {
        mesh.RecalculateNormals();
    }
}
