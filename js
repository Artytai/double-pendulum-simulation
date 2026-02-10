/* ========= math ========= */

const M = {
  rad: d => d * Math.PI / 180,
  deg: r => r * 180 / Math.PI,
  sgn: v => (v < 0 ? -1 : v > 0 ? 1 : 0),
  clamp: (x, a, b) => Math.max(a, Math.min(b, x)),
  v: (x = 0, y = 0) => ({ x, y }),
  add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y }),
  sub: (a, b) => ({ x: a.x - b.x, y: a.y - b.y }),
  mul: (a, s) => ({ x: a.x * s, y: a.y * s }),
  dot: (a, b) => a.x * b.x + a.y * b.y,
  len: a => Math.hypot(a.x, a.y),
  norm: a => { const L = Math.hypot(a.x, a.y) || 1; return { x: a.x / L, y: a.y / L }; },
  rot: (a, ang) => ({ x: a.x * Math.cos(ang) - a.y * Math.sin(ang), y: a.x * Math.sin(ang) + a.y * Math.cos(ang) }),
  mat2: (a,b,c,d)=>({a,b,c,d}),
  mulM2V: (M2, v) => ({ x: M2.a*v.x + M2.b*v.y, y: M2.c*v.x + M2.d*v.y }),
  inv2: (M2) => { const det = M2.a*M2.d - M2.b*M2.c || 1e-12; return { a: M2.d/det, b: -M2.b/det, c: -M2.c/det, d: M2.a/det }; },
  T: v => ({x:v.y, y:v.x}), // not transpose; placeholder if needed
  zeros: n => Array.from({length:n}).fill(0),
  clampAbs: (x, m) => Math.max(-m, Math.min(m, x)),
};

/* ========= deterministic PRNG (for reproducible ensemble) ========= */
function LCG(seed=1234567){
  let s = seed >>> 0;
  return ()=> (s = (1664525*s + 1013904223) >>> 0) / 0xffffffff;
}

/* ========= integrators ========= */
// All integrate a first-order system y' = f(t,y). For rigid angles y=[th1,th2,w1,w2]. For Cartesian strings, we use custom step below.
const INTEGRATORS = {
  euler(f, t, y, dt) {
    const k = f(t, y);
    return [t + dt, y.map((v, i) => v + dt * k[i])];
  },
  semi(f, t, y, dt) {
    // Treat last half as "velocities" if structure matches; otherwise similar to symplectic kick-drift
    const dy = f(t, y);
    const n = y.length;
    const m = Math.floor(n/2);
    const y2 = y.slice();
    for (let i= m; i<n; i++) y2[i] = y[i] + dt * dy[i];     // update velocities first
    const dy2 = f(t, y2);
    for (let i= 0; i<m; i++) y2[i] = y[i] + dt * dy2[i];    // then positions with new v
    return [t + dt, y2];
  },
  verlet(f, t, y, dt) {
    // If y splits [q, v], we need accelerations a(q) as last half of f.
    // We'll attempt: compute a from f(t,y) (assumes structure [dq/dt = v, dv/dt = a])
    const n = y.length, m = Math.floor(n/2);
    const dy = f(t, y);
    const q = y.slice(0, m), v = y.slice(m);
    const a = dy.slice(m);
    const qn = q.map((qi, i) => qi + v[i]*dt + 0.5*a[i]*dt*dt);
    const yTemp = qn.concat(v);         // keep old v for a2
    const dy2 = f(t + dt, yTemp);
    const a2 = dy2.slice(m);
    const vn = v.map((vi, i) => vi + 0.5*(a[i] + a2[i])*dt);
    return [t + dt, qn.concat(vn)];
  },
  rk2(f, t, y, dt) {
    const k1 = f(t, y);
    const y2 = y.map((v,i)=>v+dt*k1[i]);
    const k2 = f(t+dt, y2);
    const yn = y.map((v,i)=>v + 0.5*dt*(k1[i] + k2[i]));
    return [t+dt, yn];
  },
  rk4(f, t, y, dt) {
    const k1 = f(t, y);
    const y2 = y.map((v, i) => v + 0.5 * dt * k1[i]);
    const k2 = f(t + 0.5 * dt, y2);
    const y3 = y.map((v, i) => v + 0.5 * dt * k2[i]);
    const k3 = f(t + 0.5 * dt, y3);
    const y4 = y.map((v, i) => v + dt * k3[i]);
    const k4 = f(t + dt, y4);
    const yn = y.map((v, i) => v + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
    return [t + dt, yn];
  },
  impl_mid(f, t, y, dt) {
    // Implicit midpoint via fixed-point iterations:
    // y_{n+1} = y_n + dt * f( t + dt/2, (y_n + y_{n+1})/2 )
    // iterate: guess y*, compute RHS, update
    const maxIt = 12, eps = 1e-8;
    let yn = y.slice();
    let guess = y.slice();
    for (let it=0; it<maxIt; it++) {
      const ym = guess.map((v,i)=>0.5*(y[i]+guess[i]));
      const k = f(t + 0.5*dt, ym);
      const newGuess = y.map((v,i)=>v + dt*k[i]);
      let diff = 0;
      for (let i=0;i<y.length;i++) diff = Math.max(diff, Math.abs(newGuess[i]-guess[i]));
      guess = newGuess;
      if (diff < eps) break;
    }
    yn = guess;
    return [t+dt, yn];
  }
};

/* ========= Rigid (Lagrangian angles) ========= */

class RigidDoublePendulum {
  constructor(opts = {}) {
    this.p = {
      g: 9.81, m1: 1, m2: 1, l1: 1, l2: 1,
      b1: 0.01, b2: 0.01, qd: 0, rm1: 0, rm2: 0,
      mag: false, magK: 600, magPos: { x: 0.5, y: 0.5 },
      wind: {x: 0, y: 0}  // only used for drag approximation
    };
    Object.assign(this.p, opts.params || {});
    this.s = { t: 0, th1: Math.PI/2, th2: Math.PI/2, w1: 0, w2: 0 };
    Object.assign(this.s, opts.state || {});
    this.dt = opts.dt ?? 0.004;
    this.intName = opts.integrator || "rk4";
    this.int = INTEGRATORS[this.intName] || INTEGRATORS.rk4;
  }
  set(o = {}) { Object.assign(this.p, o); }
  setDt(dt) { this.dt = dt; }
  setIntegrator(name) { this.intName = name; this.int = INTEGRATORS[name] || INTEGRATORS.rk4; }
  setState(s) { Object.assign(this.s, s); this.s.t = 0; }

  pos() {
    const { l1, l2 } = this.p, { th1, th2 } = this.s;
    const x1 = l1 * Math.sin(th1), y1 = l1 * Math.cos(th1);
    const x2 = x1 + l2 * Math.sin(th2), y2 = y1 + l2 * Math.cos(th2);
    return { x1, y1, x2, y2 };
  }

  energy() {
    const { m1, m2, l1, l2, g } = this.p, { th1, th2, w1, w2 } = this.s;
    const d = th2 - th1;
    const T = 0.5 * m1 * l1*l1 * w1*w1 +
              0.5 * m2 * (l1*l1*w1*w1 + l2*l2*w2*w2 + 2*l1*l2*w1*w2*Math.cos(d));
    const V = (m1+m2)*g*l1*(1-Math.cos(th1)) + m2*g*l2*(1-Math.cos(th2));
    return { T, V, E: T + V };
  }

  _alphas(th1, th2, w1, w2) {
    const { g, m1, m2, l1, l2, b1, b2, qd, rm1, rm2, mag, magK, magPos, wind } = this.p;
    const D = th2 - th1, sD = Math.sin(D), cD = Math.cos(D), eps = 1e-9;
    const den = (m1 + m2) * l1 - m2 * l1 * cD * cD + eps;
    const den2 = (l2 / l1) * den + eps;

    let a1 = ( m2*l1*w1*w1*sD*cD + m2*g*Math.sin(th2)*cD + m2*l2*w2*w2*sD - (m1+m2)*g*Math.sin(th1) ) / den;
    let a2 = ( -m2*l2*w2*w2*sD*cD + (m1+m2)*g*Math.sin(th1)*cD - (m1+m2)*l1*w1*w1*sD - (m1+m2)*g*Math.sin(th2) ) / den2;

    const I1 = ((m1+m2) + rm1/3)*l1*l1 + eps;  // crude inclusion of rod mass inertia
    const I2 = (m2 + rm2/3)*l2*l2 + eps;

    a1 += -(b1 * w1) / I1;
    a2 += -(b2 * w2) / I2;

    if (qd !== 0) {
      // approximate bob linear speeds, include "wind" relative velocity
      const v1 = Math.hypot(l1*w1 - wind.x, -wind.y);
      const v2 = Math.hypot(l1*w1 + l2*w2 - wind.x, -wind.y);
      const tau1 = -qd * v1 * v1 * M.sgn(w1);
      const tau2 = -qd * v2 * v2 * M.sgn(w2);
      a1 += tau1 / I1; a2 += tau2 / I2;
    }

    if (mag) {
      const x1 = l1*Math.sin(th1), y1=l1*Math.cos(th1);
      const x2 = x1 + l2*Math.sin(th2), y2 = y1 + l2*Math.cos(th2);
      const r1 = {x: x1 - magPos.x, y: y1 - magPos.y};
      const r2 = {x: x2 - magPos.x, y: y2 - magPos.y};
      const R1 = Math.max(0.08, Math.hypot(r1.x, r1.y));
      const R2 = Math.max(0.08, Math.hypot(r2.x, r2.y));
      const F1 = -magK / (R1*R1);
      const F2 = -magK / (R2*R2);
      const t1 = F1 * (r1.x*Math.cos(th1) - r1.y*Math.sin(th1));
      const t2 = F2 * (r2.x*Math.cos(th2) - r2.y*Math.sin(th2));
      a1 += t1 / I1; a2 += t2 / I2;
    }
    return [a1, a2];
  }

  deriv(t, y) {
    const [th1, th2, w1, w2] = y;
    const [a1, a2] = this._alphas(th1, th2, w1, w2);
    return [w1, w2, a1, a2];
  }

  stepRaw(dt) {
    const y = [this.s.th1, this.s.th2, this.s.w1, this.s.w2];
    const [tn, yn] = this.int((t,y)=>this.deriv(t,y), this.s.t, y, dt);
    [this.s.t, this.s.th1, this.s.th2, this.s.w1, this.s.w2] = [tn, ...yn];
  }

  step() { this.stepRaw(this.dt); }
  stepWith(dt) { this.stepRaw(dt); }

  sample() {
    const p = this.pos(), e = this.energy();
    return { t: this.s.t, th1: this.s.th1, th2: this.s.th2, w1: this.s.w1, w2: this.s.w2, ...p, ...e };
  }
}

/* ========= Strings (Cartesian + constraints): two point masses, two constraints =========
   q = [x1 y1 x2 y2]^T, v = dq/dt, M = diag(m1, m1, m2, m2)
   g1(q) = x1^2 + y1^2 - l1^2 = 0
   g2(q) = (x2-x1)^2 + (y2-y1)^2 - l2^2 = 0
   J = ∂g/∂q is 2x4
   Equations: M a = F + Jᵀ λ
   Constraint acceleration: J a = -Ĵ v (or with stabilization terms)
   Solve for a, λ via block system:
   | M   -Jᵀ | | a | = | F |
   | J    0  | | λ |   | -Ĵv - α g - β J v |
   Use Baumgarte stabilization (α, β), optional SHAKE projection after step.
*/

class StringConstraintSystem {
  constructor(opts={}) {
    this.p = {
      g: 9.81, m1: 1, m2: 1, l1: 1.2, l2: 1.2, qd: 0.02,
      mag: false, magK: 600, magPos: {x:0.5,y:0.5}, wind:{x:0,y:0},
      collide: false, rest: 0.8, mu: 0.25, obstacle: {x:0.15, y:1.25, r:0.16}, groundY: 2.0,
      baumgarte: true, beta: 0.2, zeta: 0.7, shakeN: 3
    };
    Object.assign(this.p, opts.params||{});
    this.s = { t:0, q:[0, this.p.l1, 0, this.p.l1+this.p.l2], v:[0,0,0,0] };
    if (opts.stateAngles) this.setAngles(opts.stateAngles.th1 ?? Math.PI/2, opts.stateAngles.th2 ?? Math.PI/2);
    this.dt = opts.dt ?? 0.004;
    this.intName = opts.integrator || "rk4";  // used for debug; we use custom step anyway
  }
  set(o={}){Object.assign(this.p,o)}
  setDt(dt){this.dt=dt}
  setAngles(th1, th2) {
    const { l1, l2 } = this.p;
    const x1 = l1 * Math.sin(th1), y1 = l1 * Math.cos(th1);
    const x2 = x1 + l2 * Math.sin(th2), y2 = y1 + l2 * Math.cos(th2);
    this.s = { t:0, q:[x1,y1,x2,y2], v:[0,0,0,0] };
  }
  angles() {
    const [x1,y1,x2,y2] = this.s.q;
    const th1 = Math.atan2(x1, y1);
    const th2 = Math.atan2(x2-x1, y2-y1);
    return { th1, th2 };
  }

  _massMatrix() {
    const { m1, m2 } = this.p;
    // diag(m1, m1, m2, m2)
    return [
      [m1, 0,  0,  0],
      [0,  m1, 0,  0],
      [0,  0,  m2, 0],
      [0,  0,  0,  m2]
    ];
  }

  _externalForces(q, v) {
    const { g, m1, m2, qd, mag, magK, magPos, wind } = this.p;
    const [x1,y1,x2,y2] = q;
    const [vx1,vy1,vx2,vy2] = v;
    const F = [0,0,0,0];

    // gravity
    F[1] += g*m1;
    F[3] += g*m2;

    // quadratic drag against relative wind (per-mass)
    if (qd !== 0) {
      const rel1 = { x: vx1 - wind.x, y: vy1 - wind.y };
      const rel2 = { x: vx2 - wind.x, y: vy2 - wind.y };
      const d1 = M.mul(M.norm(rel1), -qd * M.len(rel1) * M.len(rel1));
      const d2 = M.mul(M.norm(rel2), -qd * M.len(rel2) * M.len(rel2));
      F[0] += d1.x; F[1] += d1.y;
      F[2] += d2.x; F[3] += d2.y;
    }

    // inverse-square magnet
    if (mag) {
      const r1 = { x: x1 - magPos.x, y: y1 - magPos.y };
      const r2 = { x: x2 - magPos.x, y: y2 - magPos.y };
      const R1 = Math.max(0.05, M.len(r1));
      const R2 = Math.max(0.05, M.len(r2));
      const f1 = M.mul(M.norm(r1), -magK/(R1*R1));
      const f2 = M.mul(M.norm(r2), -magK/(R2*R2));
      F[0] += f1.x; F[1] += f1.y;
      F[2] += f2.x; F[3] += f2.y;
    }

    return F;
  }

  _constraints(q) {
    const { l1, l2 } = this.p;
    const [x1,y1,x2,y2] = q;
    const g1 = x1*x1 + y1*y1 - l1*l1;
    const dx = x2-x1, dy=y2-y1;
    const g2 = dx*dx + dy*dy - l2*l2;
    return [g1, g2];
  }

  _J(q) {
    const [x1,y1,x2,y2] = q;
    // g1 = x1^2 + y1^2 - l1^2 => ∂g1/∂q = [2x1, 2y1, 0, 0]
    // g2 = (x2-x1)^2 + (y2-y1)^2 - l2^2 => ∂g2/∂q = [-2(x2-x1), -2(y2-y1), 2(x2-x1), 2(y2-y1)]
    const dx = x2-x1, dy = y2-y1;
    return [
      [2*x1, 2*y1, 0,   0  ],
      [-2*dx, -2*dy, 2*dx, 2*dy]
    ];
  }

  _Jv(q, v) {
    // time derivative dot(g) = J v
    const J = this._J(q);
    return [
      J[0][0]*v[0] + J[0][1]*v[1] + J[0][2]*v[2] + J[0][3]*v[3],
      J[1][0]*v[0] + J[1][1]*v[1] + J[1][2]*v[2] + J[1][3]*v[3]
    ];
  }

  _Jdotv(q, v) {
    // For polynomial constraints, Ĵ v = d/dt (J) v ≈ numerical or analytic:
    // g1: J1 = [2x1, 2y1, 0, 0], so Ĵ1 v = [2vx1, 2vy1, 0, 0]·v ? Not exactly; we need (dJ/dt)*v, but for linear-in-q J, (dJ/dt) depends on qdot.
    // Instead: d^2 g/dt^2 = J a + Ĵ v = second derivative of constraint. We'll compute Ĵ v analytically:
    // g1 = x1^2 + y1^2 - l1^2 -> d/dt g1 = 2x1 vx1 + 2y1 vy1, then second: 2vx1 vx1 + 2x1 ax1 + 2vy1 vy1 + 2y1 ay1
    // So Ĵ v (for g1) = 2(vx1^2 + vy1^2)
    // g2: (x2-x1)^2 + (y2-y1)^2 -> d/dt = 2(dx)(vx2-vx1) + 2(dy)(vy2-vy1)
    // second: 2(vdx^2 + vdy^2) + 2(dx)(ax2-ax1) + 2(dy)(ay2-ay1)
    // So Ĵ v term (pieces free of a) = 2[(vx2-vx1)^2 + (vy2-vy1)^2]
    const [x1,y1,x2,y2] = q; const [vx1,vy1,vx2,vy2] = v;
    const g1 = 2*(vx1*vx1 + vy1*vy1);
    const vdx = vx2 - vx1, vdy = vy2 - vy1;
    const g2 = 2*(vdx*vdx + vdy*vdy);
    return [g1, g2];
  }

  _solveBlock(Mm, J, rhsF, rhsC) {
    // Solve:
    // [ M  -J^T ] [a] = [F]
    // [ J    0  ] [λ]   [rhsC]
    // Use Schur complement: (J M^{-1} J^T) λ = J M^{-1} F - rhsC
    const Minv = invDiag4(Mm);
    const JT = transpose2x4(J);
    // S = J Minv JT (2x2)
    const S = mul2x4_4x2(J, Minv, JT);
    // b = J Minv F - rhsC
    const MinvF = mul4x4_4x1(Minv, rhsF);    // 4
    const JMinvF = mul2x4_4x1(J, MinvF);     // 2
    const b = [ JMinvF[0] - rhsC[0], JMinvF[1] - rhsC[1] ];
    const Sinv = inv2x2(S);
    const lam = mul2x2_2x1(Sinv, b);         // λ
    // a = Minv (F - JT λ)
    const JTlam = mul4x2_2x1(JT, lam);
    const rhsA = [ rhsF[0]-JTlam[0], rhsF[1]-JTlam[1], rhsF[2]-JTlam[2], rhsF[3]-JTlam[3] ];
    const a = mul4x4_4x1(Minv, rhsA);
    return { a, lam };
  }

  step() {
    const dt = this.dt;
    const { baumgarte, beta, zeta, shakeN } = this.p;
    const m = this._massMatrix();
    const F = this._externalForces(this.s.q, this.s.v); // size 4

    const g = this._constraints(this.s.q); // size 2
    const J = this._J(this.s.q);          // 2x4
    const Jv = this._Jv(this.s.q, this.s.v);
    const Jdotv = this._Jdotv(this.s.q, this.s.v);

    let rhsC = [ -Jdotv[0], -Jdotv[1] ];
    if (baumgarte) {
      // add α g + β J v with α=β*ζ? We'll use standard: rhsC = -Jdotv - 2 ζ ω J v - ω^2 g
      // map beta->ω, zeta->ζ
      const omega = beta, zt = zeta;
      rhsC[0] += - 2*zt*omega*Jv[0] - (omega*omega)*g[0];
      rhsC[1] += - 2*zt*omega*Jv[1] - (omega*omega)*g[1];
    }

    const sol = this._solveBlock(m, J, F, rhsC);
    const a = sol.a;

    // integrate (semi-implicit Euler is robust)
    const v2 = [ this.s.v[0] + dt*a[0], this.s.v[1] + dt*a[1], this.s.v[2] + dt*a[2], this.s.v[3] + dt*a[3] ];
    const q2 = [ this.s.q[0] + dt*v2[0], this.s.q[1] + dt*v2[1], this.s.q[2] + dt*v2[2], this.s.q[3] + dt*v2[3] ];

    // SHAKE-style projection to satisfy constraints
    this._shakeProject(q2, v2, shakeN);

    this.s.q = q2; this.s.v = v2; this.s.t += dt;

    if (this.p.collide) this._collisions();
  }

  stepWith(dt) { const old = this.dt; this.dt = dt; this.step(); this.dt = old; }

  _shakeProject(q, v, iters=3) {
    const { l1, l2 } = this.p;
    for (let k=0; k<iters; k++) {
      // position correction to satisfy g(q)=0
      const [x1,y1,x2,y2] = q;
      const g1 = x1*x1 + y1*y1 - l1*l1;
      if (Math.abs(g1)Below are **three complete files**—`index.html`, `style.css`, and a **1,000+ line** `app.js`—focused on **deep physics** while keeping a clean, beautiful interface.  
The JavaScript implements:

- **Full Lagrangian rigid double pendulum** (angles from vertical, y-down).
- **Constraint-aware string (rope) mode** with **tension detection** and **slack/free-flight** transitions (inequality constraints turned on/off).
- **Spring mode** (rods as linear springs with damping).
- **Generalized forces**: viscous joint damping, **quadratic air drag** projected to generalized coordinates, optional **magnetic inverse-square field**.
- **Mass matrix** with optional **distributed rod masses** (slender rods: adds \(I=\frac{1}{3}m_r \ell^2\) to diagonal inertia).
- **Integrators**: Explicit Euler, Semi-implicit (symplectic) Euler, Velocity Verlet, RK4, and **adaptive RK45 (Dormand–Prince)** with **energy guard**.
- **Tensions** computed via kinematics → Newton’s second law in Cartesian; used for rope taut/slack switching.
- **Small-angle linearization** around equilibrium with **numeric Jacobian** and **eigendecomposition** (normal modes & small-angle frequencies).
- **Lyapunov analysis** (two methods): (1) twin trajectory divergence and (2) **Benettin algorithm** integrating **variational equations** with **finite-difference Jacobian** + periodic renormalization.
- **Poincaré map** (θ₁ upward zero-cross strobe of (θ₂, ω₂)).
- **Ensemble integration** (multi‑clone small perturbations).
- **Visuals**: smooth trails, phase plot, energy plot, strobe plot; light/dark theme; clean HUD with **FPS** and **λ**.

> Open `index.html` locally. No external libs.

---

## `index.html`
```html
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Double Pendulum — Advanced Physics Lab</title>
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <style id="pre">:root{--bg:#0d1021;--fg:#e9ecff}</style>
  <link rel="stylesheet" href="style.css" />
</head>
<body>
  <header class="top">
    <div class="brand">Double Pendulum — Advanced Physics Lab</div>
    <div class="tags">Lagrangian • Rope Constraints • Springs • RK45 • Energy Guard • Poincaré • Lyapunov • Ensemble</div>
  </header>

  <main class="layout">
    <section class="left">
      <canvas id="view" width="1200" height="720"></canvas>
      <div class="overlay">
        <div class="hud">
          <span>t <b id="t">0.000</b></span>
          <span>θ₁ <b id="th1">0.000</b></span>
          <span>θ₂ <b id="th2">0.000</b></span>
          <span>ω₁ <b id="w1">0.000</b></span>
          <span>ω₂ <b id="w2">0.000</b></span>
          <span>T <b id="Te">0.000</b></span>
          <span>V <b id="Ve">0.000</b></span>
          <span>E <b id="Ee">0.000</b></span>
          <span>FPS <b id="fps">0</b></span>
          <span>λ <b id="lyap">—</b></span>
        </div>
        <div class="legend">
          <div><i class="dot a"></i> base</div>
          <div><i class="dot b"></i> twin</div>
          <div><i class="dot c"></i> ensemble</div>
          <div><i class="dot d"></i> magnet (R‑Click)</div>
        </div>
      </div>
    </section>

    <section class="right">
      <div class="card">
        <div class="row">
          <button id="play">▶</button>
          <button id="pause">⏸</button>
          <button id="step">⏭</button>
          <button id="reset">⟲</button>
          <button id="random">🎲</button>
          <button id="shot">📷 PNG</button>
          <button id="theme">🌓 Theme</button>
        </div>
        <div class="row">
          <label>Mode
            <select id="mode">
              <option value="rigid">Rigid</option>
              <option value="string">String (rope)</option>
              <option value="spring">Springs</option>
            </select>
          </label>
          <label>Integrator
            <select id="integrator">
              <option value="rk45">RK45 (adaptive)</option>
              <option value="rk4">RK4</option>
              <option value="verlet">Verlet</option>
              <option value="semi">Semi-Implicit</option>
              <option value="euler">Euler</option>
            </select>
          </label>
          <label>dt
            <input id="dt" type="number" step="0.001" value="0.004" />
          </label>
          <label>Slow
            <input id="slow" type="checkbox" />
          </label>
        </div>
      </div>

      <div class="card">
        <div class="grid2">
          <label>g <input id="g" type="range" min="0.1" max="25" step="0.01" value="9.81" /><span id="g_v">9.81</span></label>
          <label>m1 <input id="m1" type="range" min="0.05" max="6" step="0.01" value="1" /><span id="m1_v">1.00</span></label>
          <label>m2 <input id="m2" type="range" min="0.05" max="6" step="0.01" value="1" /><span id="m2_v">1.00</span></label>
          <label>l1 <input id="l1" type="range" min="0.2" max="3" step="0.01" value="1" /><span id="l1_v">1.00</span></label>
          <label>l2 <input id="l2" type="range" min="0.2" max="3" step="0.01" value="1" /><span id="l2_v">1.00</span></label>
          <label>rod m₁ <input id="rm1" type="range" min="0" max="3" step="0.01" value="0" /><span id="rm1_v">0.00</span></label>
          <label>rod m₂ <input id="rm2" type="range" min="0" max="3" step="0.01" value="0" /><span id="rm2_v">0.00</span></label>
          <label>joint b₁ <input id="b1" type="range" min="0" max="0.5" step="0.001" value="0.01" /><span id="b1_v">0.010</span></label>
          <label>joint b₂ <input id="b2" type="range" min="0" max="0.5" step="0.001" value="0.01" /><span id="b2_v">0.010</span></label>
          <label>quad drag <input id="qd" type="range" min="0" max="0.8" step="0.001" value="0" /><span id="qd_v">0.000</span></label>
          <label>px/m <input id="scale" type="range" min="80" max="360" step="1" value="160" /><span id="scale_v">160</span></label>
          <label>trail <input id="trail" type="range" min="50" max="4000" step="10" value="1800" /><span id="trail_v">1800</span></label>
        </div>
      </div>

      <div class="card">
        <div class="grid2">
          <label>θ₁ <input id="ic_th1" type="number" step="0.001" value="1.5708" /></label>
          <label>θ₂ <input id="ic_th2" type="number" step="0.001" value="1.5708" /></label>
          <label>ω₁ <input id="ic_w1" type="number" step="0.001" value="0" /></label>
          <label>ω₂ <input id="ic_w2" type="number" step="0.001" value="0" /></label>
        </div>
        <div class="row">
          <button id="applyIC">Apply</button>
          <label><input id="dragSet" type="checkbox" /> Drag‑to‑Set</label>
          <label><input id="grab2" type="checkbox" checked /> Grab Bob 2</label>
        </div>
      </div>

      <div class="card">
        <div class="row">
          <label><input id="twin" type="checkbox" /> Twin</label>
          <label>Δ° <input id="twinDelta" type="number" step="0.0001" value="0.001" /></label>
          <label><input id="ensemble" type="checkbox" /> Ensemble (8)</label>
          <label><input id="showPhase" type="checkbox" checked /> Phase</label>
          <label><input id="showEnergy" type="checkbox" checked /> Energy</label>
          <label><input id="showStrobe" type="checkbox" checked /> Strobe</label>
        </div>
        <div class="row">
          <label><input id="mag" type="checkbox" /> Magnet</label>
          <label>k <input id="magK" type="range" min="-2000" max="2000" step="10" value="600" /> <span id="magK_v">600</span></label>
          <label>spring k <input id="sk" type="range" min="0" max="800" step="1" value="260" /> <span id="sk_v">260</span></label>
          <label>spring c <input id="sc" type="range" min="0" max="10" step="0.05" value="1" /> <span id="sc_v">1.00</span></label>
        </div>
        <canvas id="phase" width="520" height="150"></canvas>
        <canvas id="energy" width="520" height="120"></canvas>
        <canvas id="strobe" width="520" height="150"></canvas>
      </div>

      <div class="card">
        <div class="row">
          <label><input id="adaptive" type="checkbox" checked /> Adaptive dt</label>
          <label>min <input id="dtMin" type="number" step="0.001" value="0.001" /></label>
          <label>max <input id="dtMax" type="number" step="0.001" value="0.02" /></label>
          <label>ΔE/E tol <input id="etol" type="number" step="0.0001" value="0.002" /></label>
        </div>
        <div class="row">
          <button id="presetClassic">Preset: Classic</button>
          <button id="presetWild">Preset: Wild</button>
          <button id="presetDamped">Preset: Damped</button>
          <button id="presetRope">Preset: Rope</button>
        </div>
        <div class="row small">Shortcuts: <code>space</code> play/pause, <code>S</code> step, <code>R</code> reset, <code>G</code> random, <code>[</code>/<code>]</code> dt, <code>P</code> PNG.</div>
      </div>
    </section>
  </main>

  <script src="app.js"></script>
</body>
</html>
