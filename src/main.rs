
mod renderer; //Used to generate output gifs
use renderer::ParticleFrame;

#[derive(Clone)]
struct Particle {

    u: f64, //Pertubation from its lattice postion
    theta: f64, //Rotation relative to the lattice
    v: f64, //Velocity
    omega: f64, //Angular velocity
    id:usize //Index in the lattice
}

struct ParticleSystem { //For simplicity mass, moment of inertia, spring constant and dipole strength are the same for all particles in a system, hence theyre defined here
    mass:f64, //Mass of each particle
    moi:f64, //Moment of inertia of each particle
    spring:f64, //spring constant for harmonic approximation, where energy is given by 0.5*k*x^2
    field:f64, //External field strength B
    field_angle:f64, //Alignment of field relative to lattice
    dipole:f64, //Dipole strength of each particle
    spacing:f64, //Lattice constant a

    rot_damping:f64, //Constants for simple first order damping
    linear_damping:f64,

    particles: Vec<Particle>,
}


impl ParticleSystem {
    fn new(mass:f64,moi:f64,spring:f64,field:f64,field_angle:f64,dipole:f64,spacing:f64, rot_damping:f64,linear_damping:f64, mut particles: Vec<Particle>) -> Self {
        particles.sort_by_key(|p| p.id);
        Self { mass, moi, spring, field, field_angle, dipole, spacing,rot_damping, linear_damping, particles }
    }

    fn get_neighbours(&self, id: usize) -> [&Particle; 2] { //Returns adjacent particles in a system, periodically
        let n = self.particles.len();

        let left = (id + n - 1) % n;
        let right = (id + 1) % n;

        [
            &self.particles[left],
            &self.particles[right],
        ]
    }

    fn get_linear_contributions(&self, id:usize) -> f64{
        let mut F=0.0;
        let q = &self.particles[id];
        for p in &self.particles{
            if p.id==id{
                continue;
            }
            let magnetic_cont=0.5*3.0*(self.spacing*(((id as f64)-(p.id as f64)))+q.u-p.u).signum()*(self.spacing*(((id as f64)-(p.id as f64)))+q.u-p.u).powf(-4.0)*((q.theta-p.theta).cos()-3.0*q.theta.cos()*p.theta.cos());

            F+=magnetic_cont;
        }
        return F;
    }
    fn get_angular_contributions(&self, id:usize) -> f64{
        let mut F=0.0;
        let q = &self.particles[id];
        for p in &self.particles{
            if p.id==id{
                continue;
            }
            let magnetic_cont=-0.5*((self.spacing*(((id as f64)-(p.id as f64)))+q.u-p.u).abs().powf(-3.0)*(-((q.theta-p.theta).sin())+3.0*(p.theta.cos())*(q.theta.sin())));
            F+=magnetic_cont;
        }
        return F;
    }

    fn step(&self,timestep:f64) -> Self { //Retuns an updated system one timestep in the future
        let mut updated = self.particles.clone();

        for p in &self.particles {
            let acc = find_acceleration(self, p.id);
            updated[p.id].v+=timestep*acc[0];
            updated[p.id].omega+=timestep*acc[1];
            updated[p.id].u+=timestep*updated[p.id].v;
            updated[p.id].theta+=timestep*updated[p.id].omega;
        }

        return ParticleSystem::new(self.mass,self.moi,self.spring,self.field,self.field_angle,self.dipole,self.spacing, self.rot_damping, self.linear_damping, updated);
    }

    fn populate(&mut self, count:usize){ //Fills an empty system with some number of default particles
        let mut ps:Vec<Particle> = vec![];
        for i in 0..count{
            let particle = Particle {
                id: i,
                u: 0.0,
                theta: 0.0,
                v: 0.0,
                omega: 0.0,
            };
            ps.push(particle);
        
        }
        self.particles=ps;
    }

}

fn find_acceleration(system: &ParticleSystem, current_particle: usize)->[f64;2]{ //Calculates the linear and rotational acceleration from the EOM
    let pos_n = system.particles[current_particle].u;
    let rot_n = system.particles[current_particle].theta;
    let v = system.particles[current_particle].v;
    let omega = system.particles[current_particle].omega;

    let neighbours = system.get_neighbours(current_particle);
    let pos_prev = neighbours[0].u;
    let rot_prev = neighbours[0].theta;
    let pos_next = neighbours[1].u;
    let rot_next = neighbours[1].theta;

    //EOM for p from Hamiltonian
    //let lin_momentum = -v*system.linear_damping-system.spring*(2.0*pos_n-pos_prev-pos_next)+3.0*system.dipole.powf(2.0)*(((rot_n-rot_prev).cos()-3.0*rot_n.cos()*rot_prev.cos())/((system.spacing + pos_n-pos_prev).powf(4.0))-((rot_next-rot_n).cos()-3.0*rot_n.cos()*rot_next.cos())/((system.spacing + pos_next-pos_n).powf(4.0)));
    let lin_momentum = -v*system.linear_damping-system.spring*(2.0*pos_n-pos_prev-pos_next)+system.get_linear_contributions(current_particle);
    //EOM for L from Hamiltonian
    //let rot_momentum = -omega*system.rot_damping-system.dipole*system.field*(rot_n-system.field_angle).sin()-3.0*system.dipole.powf(2.0)*(((rot_next-rot_n).sin()+3.0*rot_n.sin()*rot_next.cos())/((system.spacing+pos_next-pos_n).powf(3.0))-((rot_n-rot_prev).sin()-3.0*rot_n.sin()*rot_prev.cos())/((system.spacing+pos_n-pos_prev).powf(3.0)));
    let rot_momentum = -omega*system.rot_damping-system.dipole*system.field*(rot_n-system.field_angle).sin()+system.get_angular_contributions(current_particle);
    return [lin_momentum/system.mass,rot_momentum/system.moi];

}

fn main()  -> Result<(), Box<dyn std::error::Error>>{
    //Currently doing everything manually...
    let mut new_system = ParticleSystem::new(
        1.0,
        1.0,
        5.0,
        10.0,
        1.0,
        0.1,
        100.0,
        0.0,
        0.0,
        vec![]);
    new_system.populate(8);
    for p in &mut new_system.particles{
        p.u=((4.0* (3.14159265358979) *(p.id as f64))/8.0).cos()*0.1
    }

    let mut history: Vec<ParticleFrame> = Vec::new();
    let t=0.01;
    let length = 8.0;
    for _ in 0..((length/t) as usize) {
        new_system=new_system.step(t);

        history.push(ParticleFrame {
            positions: new_system.particles.iter()
                .map(|p| p.u)
                .collect(),

            angles: new_system.particles.iter()
                .map(|p| p.theta)
                .collect(),
});
    }
    renderer::render(&history, "particles.gif",t, 10)?;
    Ok(())
}
