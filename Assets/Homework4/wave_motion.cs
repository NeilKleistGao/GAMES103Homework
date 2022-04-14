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

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;

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
				new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping;
				
				for (int k = 0; k < 4; ++k) {
					int ni = i + next[k].x, nj = j + next[k].y;
					if (ni > -1 && nj > -1 && ni < size && nj < size) {
						new_h[i, j] += rate * (h[ni, nj] - h[i, j]);
					}
				}
			}
		}


		// int li = size - 1, ui = 0, lj = size - 1, uj = 0; 
		// Step 2: Block->Water coupling
		// calculate low_h.
		foreach (GameObject cube in cubes) {
			Vector3 cubePos = cube.transform.position;
			int li = Mathf.Max(0, Mathf.FloorToInt((cubePos.x + size * 0.05f) * 10.0f) - 3);
			int ui = Mathf.Min(size - 1, Mathf.CeilToInt((cubePos.x + size * 0.05f) * 10.0f) + 3);
			int lj = Mathf.Max(0, Mathf.FloorToInt((cubePos.y + size * 0.05f) * 10.0f) - 3);
			int uj = Mathf.Min(size - 1, Mathf.CeilToInt((cubePos.y + size * 0.05f) * 10.0f) + 3);

			// li = Mathf.Min(li, sli); ui = Mathf.Max(ui, sui);
			// lj = Mathf.Min(lj, slj); uj = Mathf.Max(uj, suj);

			for (int i = li; i <= ui; ++i) {
				for (int j = lj; j <= uj; ++j) {
					low_h[i, j] = new_h[i, j];
					Vector3 pos = new Vector3(i*0.1f-size*0.05f, 0, j*0.1f-size*0.05f);

					Vector3 direction = (cube.transform.position - pos);
					Ray ray = new Ray(pos, direction.normalized);
					Vector3 temp = direction;
					temp.y = 0;
					float cosTheta = Vector3.Dot(direction.normalized, temp.normalized);
					float maxDistance = 0.5f / cosTheta;

					RaycastHit[] hits = Physics.RaycastAll(ray, maxDistance);
					foreach (RaycastHit hit in hits) {
						low_h[i, j] = -hit.distance / maxDistance;
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

		//Step 4: Water->Block coupling.
		//More TODO here.
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
	}
}
