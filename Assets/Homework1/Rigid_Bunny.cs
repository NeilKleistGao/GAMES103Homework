using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;					// for collision

	Mesh mesh = null;
	Vector3[] vertices = null;


	// Use this for initialization
	void Start () 
	{		
		mesh = GetComponent<MeshFilter>().mesh;
		vertices = mesh.vertices;

		float m = 1.0f;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

	// In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		Vector3 collisionVertex = Vector3.zero;
		int collisionCount = 0;
		Matrix4x4 R = GetRotationMatrix(transform.rotation);
		Matrix4x4 model = transform.localToWorldMatrix;

		foreach (Vector3 vet in vertices) {
			Vector4 h = new Vector4(vet.x, vet.y, vet.z, 1);
			Vector3 wv = model * h;
			Vector3 rad = R * h;
			Vector3 velocity = v + Vector3.Cross(w, rad);

			float SDF = Vector3.Dot(wv - P, N.normalized);
			if (SDF < 0 && Vector3.Dot(velocity, N) < 0) {
				collisionVertex += vet;
				++collisionCount;
			}
		}

		if (collisionCount > 0) {
			collisionVertex /= collisionCount;
			Vector4 h = new Vector4(collisionVertex.x, collisionVertex.y, collisionVertex.z, 1);
			Vector3 rad = R * h;
			Vector3 velocity = v + Vector3.Cross(w, rad);

			Vector3 vn = Vector3.Dot(velocity, N.normalized) * N.normalized;
			Vector3 vt = v - vn;

			Vector3 nextVN = -restitution * vn;
			float alpha = Mathf.Max(1 - 0.9f * (1 + restitution) * vn.magnitude / vt.magnitude, 0);
			Vector3 nextVT = alpha * vt;
			Vector3 nextV = nextVN + nextVT;

			Matrix4x4 I = R * I_ref * R.transpose;
			Matrix4x4 K = Matrix4x4.identity;
			K[0, 0] /= mass; K[1, 1] /= mass; K[2, 2] /= mass;
			Matrix4x4 temp = Get_Cross_Matrix(rad) * I.inverse * Get_Cross_Matrix(rad);
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					K[i, j] -= temp[i, j];
				}
			}

			Vector4 deltaV = nextV - velocity;
			deltaV.w = 0;
			Vector3 J = K.inverse * deltaV;

			v += (1.0f / mass) * J;
			Vector4 deltaW = I.inverse * Vector3.Cross(rad, J);
			w += new Vector3(deltaW.x, deltaW.y, deltaW.z);
		}
	}

	Matrix4x4 GetRotationMatrix(Quaternion q) {
		Matrix4x4 R = Matrix4x4.identity;
		R[0, 0] = q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z;
		R[0, 1] = 2 * (q.x * q.y - q.w * q.z);
		R[0, 2] = 2 * (q.x * q.z + q.w * q.y);
		R[1, 1] = q.w * q.w - q.x * q.x + q.y * q.y - q.z * q.z;
		R[1, 0] = 2 * (q.x * q.y + q.w * q.z);
		R[1, 2] = 2 * (q.y * q.z - q.w * q.x);
		R[2, 2] = q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z;
		R[2, 0] = 2 * (q.x * q.z - q.w * q.y);
		R[2, 1] = 2 * (q.y * q.z + q.w * q.x);

		return R;
	}

	// Update is called once per frame
	void Update () 
	{
		//Game Control
		if(launched && Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
			restitution = 0.5f;
			launched=false;
		}
		if(!launched && Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
			launched=true;
		}

		if (!launched) {
			return;
		}

		// Part I: Update velocities
		Vector3 force = new Vector3(0, -9.8f, 0);
		v = v + dt * force;
		v *= linear_decay;
		w *= angular_decay;

		// w = new Vector3(0, 100, 0); // for test
		// No need for tau update since the unique force is gravity.

		// Part II: Collision Impulse
		Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		//Update linear status
		Vector3 x    = transform.position;
		x = x + dt * v;

		//Update angular status
		Quaternion q = transform.rotation;
		Quaternion deltaQ = new Quaternion();
		deltaQ.w = 0;
		deltaQ.x = dt * 0.5f * w.x;
		deltaQ.y = dt * 0.5f * w.y;
		deltaQ.z = dt * 0.5f * w.z;
		deltaQ = deltaQ * q;
		q.x += deltaQ.x; q.y += deltaQ.y; q.z += deltaQ.z; q.w += deltaQ.w;

		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
