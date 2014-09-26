#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

#ifndef INFINITY
#define INFINITY 1e8
#endif

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#define MAX_RAY_DEPTH 5
#define BIAS 1e-4

typedef float FloatType;

class Vector3
{
public:
	FloatType x, y, z;
	Vector3() : x(FloatType(0)), y(FloatType(0)), z(FloatType(0)) {}
	Vector3(FloatType xx) : x(xx), y(xx), z(xx) {}
	Vector3(FloatType xx, FloatType yy, FloatType zz) : x(xx), y(yy), z(zz) {}
	Vector3& normalize()
	{
		FloatType nor2 = lengthSqr();
		if (nor2 > 0) {
			FloatType invNor = 1 / sqrt(nor2);
			x *= invNor, y *= invNor, z *= invNor;
		}
		return *this;
	}
	Vector3 operator * (const FloatType &f) const { return Vector3(x * f, y * f, z * f); }
	Vector3 operator * (const Vector3 &v) const { return Vector3(x * v.x, y * v.y, z * v.z); }
	FloatType dot(const Vector3 &v) const { return x * v.x + y * v.y + z * v.z; }
	Vector3 operator - (const Vector3 &v) const { return Vector3(x - v.x, y - v.y, z - v.z); }
	Vector3 operator + (const Vector3 &v) const { return Vector3(x + v.x, y + v.y, z + v.z); }
	Vector3& operator += (const Vector3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
	Vector3& operator *= (const Vector3 &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	Vector3 operator - () const { return Vector3(-x, -y, -z); }
	FloatType lengthSqr() const { return x * x + y * y + z * z; }
	FloatType length() const { return sqrt(lengthSqr()); }
};

class Intersecter
{
public:
   virtual bool Intersect(const Vector3 &rayOrig, const Vector3 &rayDir, FloatType &t1, FloatType &f2) const = 0;
   virtual Vector3 Normal(const Vector3 &intersectPoint) const = 0;
   virtual FloatType GetTransparency() const = 0;
   virtual FloatType GetReflection() const = 0;
   virtual const Vector3 &GetSurfaceColor() const = 0;
   virtual const Vector3 &GetEmissionColor() const = 0;
};


FloatType mix(const FloatType &a, const FloatType &b, const FloatType &mix)
{
	return b * mix + a * (FloatType(1) - mix);
}

class Sphere : public Intersecter
{
public:
   Vector3 center;                         /// position of the sphere
   FloatType radius, radius2;                      /// sphere radius and radius^2
   Vector3 surfaceColor, emissionColor;    /// surface color and emission (light)
   FloatType transparency, reflection;             /// surface transparency and reflectivity
   Sphere(const Vector3 &c, const FloatType &r, const Vector3 &sc, 
          const FloatType &refl = 0, const FloatType &transp = 0, const Vector3 &ec = 0) : 
      center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
      transparency(transp), reflection(refl)
   {}
   // compute a ray-sphere intersection using the geometric solution
   virtual bool Intersect(const Vector3 &rayorig, const Vector3 &raydir, FloatType &t0, FloatType &t1) const
   {
      Vector3 l = center - rayorig;
      FloatType tca = l.dot(raydir);
      if (tca < 0) return false;
      FloatType d2 = l.dot(l) - tca * tca;
      if (d2 > radius2) return false;
      FloatType thc = sqrt(radius2 - d2);
      t0 = tca - thc;
      t1 = tca + thc;
      
      return true;
   }

   virtual Vector3 Normal(const Vector3 &intersectPoint) const
   {
      return intersectPoint - this->center;
   }
   
   virtual FloatType GetTransparency() const { return transparency; }
   virtual FloatType GetReflection() const { return reflection; }
   virtual const Vector3 &GetSurfaceColor() const { return surfaceColor; }
   virtual const Vector3 &GetEmissionColor() const { return emissionColor; }
};

Vector3 trace(const Vector3 &rayOrig, const Vector3 &rayDir, 
              const std::vector<Intersecter *> &objects, const int &depth);

Vector3 getRefraction(const Vector3 &rayDir, const Vector3 &phit, const Vector3 &nhit, bool inside, const std::vector<Intersecter *> &objects, int depth)
{
   Vector3 refraction = 0;
   FloatType ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
   FloatType cos = -nhit.dot(rayDir);
   FloatType k = 1 - eta * eta * (1 - cos * cos);
   Vector3 refractionDir = rayDir * eta + nhit * (eta *  cos - sqrt(k));
   refractionDir.normalize();
   refraction = trace(phit - nhit * BIAS, refractionDir, objects, depth);

   return refraction;
}

Vector3 getReflection(const Vector3 &rayDir, const Vector3 &phit, const Vector3 &nhit, const std::vector<Intersecter *> &objects, int depth, FloatType &fresnelEffect)
{
   FloatType facingratio = -rayDir.dot(nhit);
   // change the mix value to tweak the effect
   fresnelEffect = mix(pow(1 - facingratio, 3), 1, 0.1); 
   // compute reflection direction (not need to normalize because all vectors
   // are already normalized)
   Vector3 refldir = rayDir - nhit * 2 * rayDir.dot(nhit);
   refldir.normalize();
   Vector3 reflection = trace(phit + nhit * BIAS, refldir, objects, depth);
   return reflection;
}

Intersecter *findClosest(const Vector3 &rayOrig, const Vector3 &rayDir, const std::vector<Intersecter *> &objects, FloatType &tnear)
{
   Intersecter *closest = NULL;
   
   // find intersection of this ray with the closest in the scene
   for (unsigned i = 0; i < objects.size(); ++i) {
      FloatType t0 = INFINITY, t1 = INFINITY;
      if (objects[i]->Intersect(rayOrig, rayDir, t0, t1)) {
         if (t0 < 0) t0 = t1;
         if (t0 < tnear) {
            tnear = t0;
            closest = objects[i];
         }
      }
   }

   return closest;
}

Vector3 trace(const Vector3 &rayOrig, const Vector3 &rayDir, 
	const std::vector<Intersecter *> &objects, const int &depth)
{
   FloatType tnear = INFINITY;
   const Intersecter *closest = findClosest(rayOrig, rayDir, objects, tnear);

   if (!closest) return Vector3(2);

   Vector3 surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
   Vector3 phit = rayOrig + rayDir * tnear; // point of intersection
   Vector3 nhit = closest->Normal(phit);
   nhit.normalize();
   
   bool inside = rayDir.dot(nhit) > 0;
   if (inside)
   {
      nhit = -nhit;
   }
   
   if ((closest->GetTransparency() > 0 || closest->GetReflection() > 0) && depth < MAX_RAY_DEPTH)
   {
      FloatType fresnelEffect;
      Vector3 reflection = getReflection(rayDir, phit, nhit, objects, depth + 1, fresnelEffect);
      
      Vector3 refraction = 0;
      // if the closest is also transparent compute refraction ray (transmission)
      if (closest->GetTransparency()) {
         refraction = getRefraction(rayDir, phit, nhit, inside, objects, depth + 1);
      }

      surfaceColor = (reflection * fresnelEffect + 
                      refraction * (1 - fresnelEffect) * closest->GetTransparency()) * closest->GetSurfaceColor();
   }
   else {
      for (unsigned i = 0; i < objects.size(); ++i) {
         if (objects[i]->GetEmissionColor().x > 0) {
            // this is a light
            Vector3 transmission = 1;
            Vector3 lightDirection = -(objects[i]->Normal(phit));
            lightDirection.normalize();
            for (unsigned j = 0; j < objects.size(); ++j) {
               if (i != j) {
                  FloatType t0, t1;
                  if (objects[j]->Intersect(phit + nhit * BIAS, lightDirection, t0, t1)) {
                     transmission = 0;
                     break;
                  }
               }
            }
            surfaceColor += closest->GetSurfaceColor() * transmission * 
               std::max(FloatType(0), nhit.dot(lightDirection)) * objects[i]->GetEmissionColor();
         }
      }
   }

   return surfaceColor + closest->GetEmissionColor();
}
