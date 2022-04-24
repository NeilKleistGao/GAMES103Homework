using UnityEngine;
using System.Collections;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.004f;
	float damping 	= 0.996f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	Vector3[] 	cube_v = new Vector3[2] { Vector3.zero, Vector3.zero };
	Vector3[] 	cube_w = new Vector3[2] { Vector3.zero, Vector3.zero };

	System.Random random = new System.Random();

	[SerializeField] private GameObject[] cubes = null;

	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}

	void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
	{		
		Vector2Int[] next = new Vector2Int[4] {
			new Vector2Int(1, 0), new Vector2Int(0, 1), new Vector2Int(-1, 0), new Vector2Int(0, -1)
		};
		// Step 1:
		// Compute new_h based on the shallow wave model.
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				low_h[i, j] = h[i, j];
				new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping;
				
				for (int k = 0; k < 4; ++k) {
					int ni = i + next[k].x, nj = j + next[k].y;
					if (ni > -1 && nj > -1 && ni < size && nj < size) {
						new_h[i, j] += rate * (h[ni, nj] - h[i, j]);
					}
				}
			}
		}


		// Step 2: Block->Water coupling
		// calculate low_h.
		int cubeIndex = 0;
		foreach (GameObject cube in cubes) {
			Vector3 cubePos = cube.transform.position;
			if (cubePos.y >= 0.5f) {
				continue;
			}

			int li = Mathf.Max(0, Mathf.FloorToInt((cubePos.x + size * 0.05f) * 10.0f) - 3);
			int ui = Mathf.Min(size - 1, Mathf.CeilToInt((cubePos.x + size * 0.05f) * 10.0f) + 3);
			int lj = Mathf.Max(0, Mathf.FloorToInt((cubePos.z + size * 0.05f) * 10.0f) - 3);
			int uj = Mathf.Min(size - 1, Mathf.CeilToInt((cubePos.z + size * 0.05f) * 10.0f) + 3);

			for (int i = li; i <= ui; ++i) {
				for (int j = lj; j <= uj; ++j) {
					Vector3 pos = new Vector3(i*0.1f-size*0.05f, -1, j*0.1f-size*0.05f);

					Vector3 direction = (cube.transform.position - pos);
					Ray ray = new Ray(pos, direction.normalized);

					RaycastHit[] hits = Physics.RaycastAll(ray, 1);
					foreach (RaycastHit hit in hits) {
						low_h[i, j] = Mathf.Max(-1.0f/hit.distance, -5);
					}
				}
			}

			// then set up b and cg_mask for conjugate gradient.
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					if (low_h[i, j] < h[i, j]) {
						cg_mask[i, j] = true;
						b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
					}
					else {
						cg_mask[i, j] = false;
						b[i, j] = 0;
						vh[i, j] = 0;
					}
				}
			}
			// Solve the Poisson equation to obtain vh (virtual height).
			// Diminish vh.
			Conjugate_Gradient(cg_mask, b, vh, li, ui, lj, uj);
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					vh[i, j] *= gamma;
				}
			}
			// Update new_h by vh.
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					for (int k = 0; k < 4; ++k) {
						int ni = i + next[k].x, nj = j + next[k].y;
						if (ni > -1 && nj > -1 && ni < size && nj < size) {
							new_h[i, j] += (vh[ni, nj] - vh[i, j]) * rate;
						}
					}
				}
			}
		}

		// Step 3
		// old_h <- h; h <- new_h;
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				old_h[i, j] = h[i, j];
				h[i, j] = new_h[i, j];
			}
		}
	}
	

	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

		// Load X.y into h.
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				h[i, j] = X[i * size + j].y;
			}
		}

		if (Input.GetKeyDown ("r")) 
		{
			// Add random water.
			int x = random.Next(0, size), y = random.Next(0, size);
			float r = Mathf.Clamp01((float)random.NextDouble()) / 2;
			int surroundCount = 0;
			if (x > 0) { 
				++surroundCount;
			}
			if (x < size - 1) { 
				++surroundCount;
			}
			if (y > 0) {
				++surroundCount;
			}
			if (y < size - 1) {
				++surroundCount;
			}

			h[x, y] += r;
			if (x > 0) {
				h[x - 1, y] -= r / surroundCount;
			}
			if (x < size - 1) {
				h[x + 1, y] -= r / surroundCount;	
			}
			if (y > 0) {
				h[x, y - 1] -= r / surroundCount;
			}
			if (y < size - 1) {
				h[x, y + 1] -= r / surroundCount;
			}
		}
	
		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
		}

		// Store h back into X.y and recalculate normal.
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				X[i * size + j].y = h[i, j];
			}
		}
		
		mesh.vertices = X;
		mesh.RecalculateNormals();

		{
			GameObject cube = cubes[1];
			Vector3 cubePos = cube.transform.position;
			int li = Mathf.Max(0, Mathf.FloorToInt((cubePos.x + size * 0.05f) * 10.0f) - 3);
			int ui = Mathf.Min(size - 1, Mathf.CeilToInt((cubePos.x + size * 0.05f) * 10.0f) + 3);
			int lj = Mathf.Max(0, Mathf.FloorToInt((cubePos.z + size * 0.05f) * 10.0f) - 3);
			int uj = Mathf.Min(size - 1, Mathf.CeilToInt((cubePos.z + size * 0.05f) * 10.0f) + 3);

			Vector3 force = new Vector3(0, -9.8f, 0);
			Vector3 torque = Vector3.zero;
			if (cubePos.y < 0.5f) {
				for (int i = li; i <= ui; ++i) {
					for (int j = lj; j <= uj; ++j) {
						if (vh[i,j] > 0) {
							Vector3 f = new Vector3(0, 9.8f * vh[i,j] * 0.02f, 0);
							force += f;

							Vector3 pos = new Vector3(i*0.1f-size*0.05f, -1, j*0.1f-size*0.05f);

							Vector3 direction = (cube.transform.position - pos);
							Ray ray = new Ray(pos, direction.normalized);
							RaycastHit[] hits = Physics.RaycastAll(ray, 1);
							foreach (RaycastHit hit in hits) {
								if (hit.transform == cube.transform) {
									torque += Vector3.Cross(hit.point - cube.transform.position, f);
								}
							}
						}
					}
				}
			}

			cube_v[1] += force * 0.004f;
			cube_v[1] *= damping;
			Vector3 res = cube.transform.position + cube_v[1] * 0.004f;
			cube.transform.position = res;

			cube_w[1] += torque * 0.004f * 1000;
			cube_w[1] *= damping;
			Quaternion q = cube.transform.rotation;
			Quaternion deltaQ = new Quaternion();
			deltaQ.w = 0;
			deltaQ.x = 0.004f * 0.5f * cube_w[1].x;
			deltaQ.y = 0.004f * 0.5f * cube_w[1].y;
			deltaQ.z = 0.004f * 0.5f * cube_w[1].z;
			deltaQ = deltaQ * q;
			q.x += deltaQ.x; q.y += deltaQ.y; q.z += deltaQ.z; q.w += deltaQ.w;
			cube.transform.rotation = q;
		}
	}
}
