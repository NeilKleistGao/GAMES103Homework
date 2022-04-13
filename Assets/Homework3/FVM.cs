using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class FVM : MonoBehaviour
{
	float dt 			= 0.003f;
    float mass 			= 1;
	float stiffness_0	= 20000.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;

	int[] 		Tet;
	int tet_number;			//The number of tetrahedra

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;
	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;

	//For Laplacian smoothing.
	Vector3[]   V_sum;
	int[]		V_num;

	SVD svd = new SVD();

    // Start is called before the first frame update
    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/Homework3/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);
        	Tet = new int[tet_number*4];

    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/Homework3/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}
        /*tet_number=1;
        Tet = new int[tet_number*4];
        Tet[0]=0;
        Tet[1]=1;
        Tet[2]=2;
        Tet[3]=3;

        number=4;
        X = new Vector3[number];
        V = new Vector3[number];
        Force = new Vector3[number];
        X[0]= new Vector3(0, 0, 0);
        X[1]= new Vector3(1, 0, 0);
        X[2]= new Vector3(0, 1, 0);
        X[3]= new Vector3(0, 0, 1);*/


        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

		for (int i = 0; i < number; ++i) {
			Force[i] = V[i] = Vector3.zero;
		}

		inv_Dm = new Matrix4x4[tet_number];
		for (int i = 0; i < tet_number; ++i) {
			inv_Dm[i] = Build_Edge_Matrix(i).inverse;
		}
    }

    Matrix4x4 Build_Edge_Matrix(int tet)
    {
    	Matrix4x4 ret=Matrix4x4.identity;
    	Vector3 v0 = X[Tet[tet * 4]], v1 = X[Tet[tet * 4 + 1]], v2 = X[Tet[tet * 4 + 2]], v3 = X[Tet[tet * 4 + 3]];
		Vector3 x10 = v1 - v0, x20 = v2 - v0, x30 = v3 - v0;

		ret[0, 0] = x10.x; ret[1, 0] = x10.y; ret[2, 0] = x10.z;
		ret[0, 1] = x20.x; ret[1, 1] = x20.y; ret[2, 1] = x20.z;
		ret[0, 2] = x30.x; ret[1, 2] = x30.y; ret[2, 2] = x30.z;

		return ret;
    }

    void _Update()
    {
    	// Jump up.
		if(Input.GetKeyDown(KeyCode.Space))
    	{
    		for(int i=0; i<number; i++)
    			V[i].y+=0.2f;
    	}

    	for(int i=0; i<number; i++)
    	{
    		Force[i] = new Vector3(0, -9.8f, 0);
    	}

    	for(int tet=0; tet<tet_number; tet++)
    	{
    		Matrix4x4 F = Build_Edge_Matrix(tet) * inv_Dm[tet];
    		Matrix4x4 G = F.transpose * F;
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					G[i, j] -= (i == j) ? 1 : 0;
					G[i, j] *= 0.5f;
				}
			}

			float trace = 0.0f;
			for (int i = 0; i < 3; ++i) {
				trace += G[i, i];
			}

    		Matrix4x4 S = Matrix4x4.zero;
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					S[i, j] = stiffness_1 * 2 * G[i, j];
					if (i == j) {
						S[i, j] += stiffness_0 * trace;
					}
				}
			}

			Matrix4x4 res = F * S * inv_Dm[tet].transpose;
			float cof = -1.0f / (6.0f * inv_Dm[tet].determinant);

			Force[Tet[tet * 4]].x -= cof * (res[0, 0] + res[0, 1] + res[0, 2]);
			Force[Tet[tet * 4]].y -= cof * (res[1, 0] + res[1, 1] + res[1, 2]);
			Force[Tet[tet * 4]].z -= cof * (res[2, 0] + res[2, 1] + res[2, 2]);
			Force[Tet[tet * 4 + 1]].x += cof * res[0, 0];
			Force[Tet[tet * 4 + 1]].y += cof * res[1, 0];
			Force[Tet[tet * 4 + 1]].z += cof * res[2, 0];
			Force[Tet[tet * 4 + 2]].x += cof * res[0, 1];
			Force[Tet[tet * 4 + 2]].y += cof * res[1, 1];
			Force[Tet[tet * 4 + 2]].z += cof * res[2, 1];
			Force[Tet[tet * 4 + 3]].x += cof * res[0, 2];
			Force[Tet[tet * 4 + 3]].y += cof * res[1, 2];
			Force[Tet[tet * 4 + 3]].z += cof * res[2, 2];
			// Force[Tet[tet * 4]] = -(Force[Tet[tet * 4 + 1]] + Force[Tet[tet * 4 + 2]] + Force[Tet[tet * 4 + 3]]);
    	}

		
    	for(int i=0; i<number; i++)
    	{
			V[i] += Force[i] * dt;
			V[i] *= damp;
    		X[i] += V[i] * dt;

			
    		if (X[i].y < -3) {
				Vector3 res = new Vector3(X[i].x, -3, X[i].z);
				V[i] += (res - X[i]) / dt;
				V[i].x = V[i].z = 0;
				X[i] = res;
			}
    	}
    }

    // Update is called once per frame
    void Update()
    {
    	for(int l=0; l<10; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals ();
    }
}
