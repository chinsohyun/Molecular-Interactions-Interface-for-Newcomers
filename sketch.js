const dt = 0.002;
const fps = 1/0.002;
// const sigma = 3.4 * Math.pow(10, -10);
// const mass = 6.63 * Math.pow(10, -26);
// const epsilon = 1.66 * Math.pow(10, -21);
const tKelvin = 100;
const tInit = 1; //1.0645968192771087, (scipy.constants.k * tKelvin) / epsilon
const rcut = 2.5;
const boxSize = 6.8/2;

let fileLines; 
let fileInput
let atomCount;
let pos;
let forces;
let velocities
let potentialE;
let isDrawing = false; 
let initialPos

let kE = [];
let pE = [];
let ha = [];

let myColor = [];

function preload() {
  importFile('https://files.cargocollective.com/c1052065/8.txt');

}

function resetSketch() {
  importFile('https://files.cargocollective.com/c1052065/8.txt');

  myColor = [];
  for (let atom = 0; atom < atomCount; atom++) { 
    myColor.push(color(random(255), random(255), random(255)));    
  }
  
}

function setup() {
  createCanvas(windowWidth, windowHeight);
  background(255);
  frameRate(fps) //FPS
  
  fileInput = createFileInput(handleFile);
  fileInput.position(5, 690);
  
  velocities = initVel();
  [ forces, potentialE ] = calculateForces(pos);
  
  for (let atom = 0; atom < atomCount; atom++) { 
    myColor.push(color(random(100,255), random(100,255), random(100,255)));    
  }
  
  initialPos = pos;
  
  
  
  //button

  let resetButton = createButton('Reset');
  resetButton.position(360, 355); // Set the position of the button
  resetButton.mousePressed(resetSketch); // Call resetSketch function on button press  
  
  let startButton = createButton('Start Drawing');
  startButton.position(360, 355+25);
  startButton.mousePressed(startDrawing);
  
  // Create a button to stop drawing
  let stopButton = createButton('Stop Drawing');
  stopButton.position(360+100, 355 +25);
  stopButton.mousePressed(stopDrawing);
}


function draw() {
    background(255);
    stroke(255)
    fill(0, 0, 150);
    rect(0,0, 350, 340);
    fill(150, 0, 0);
    rect(0,340, 350, 340);
    stroke(0)
    fill(255)
    
    textFont();
    textSize(14);
    textAlign(LEFT);
    text("Top view", 10, 20);
    text("Front view", 10, 360);

  
    const posOut = new Array(atomCount);
    const velHalf = new Array(atomCount);
    const velOut = new Array(atomCount);
  
    for (let atom = 0; atom < atomCount; atom++) { 
        const pInit = pos[atom];
        const vInit = velocities[atom];
        const fInit = forces[atom];
        const vHalf = updateVel(vInit, fInit);
        const pNew = updatePos(pInit, vHalf);
        velHalf[atom] = vHalf;
        posOut[atom] = pNew;
    }

    const [fOut, pot] = calculateForces(posOut);

    for (let atom = 0; atom < atomCount; atom++) {
        const velocity = updateVel(velHalf[atom], fOut[atom]);
        velOut[atom] = velocity
    }
       
    pos = posOut
    velocities = velOut
    forces = fOut
    potentialE = pot

    let newVelArray = extractVelocities(velocities);
    let newVelMean = calculateMean(newVelArray);
    let newVelStd = calculateStandardDeviation(newVelArray, newVelMean);
    let kineticE = (3 / 2) * atomCount * newVelStd ** 2;
  
    let hamilton = kineticE + potentialE

    kE.push(kineticE);
    pE.push(potentialE);
    ha.push(hamilton);  
    if (kE.length < windowWidth - 450) {
    kE.push(kineticE);
    pE.push(potentialE);
    ha.push(hamilton);

  } else {
    kE = [];
    pE = [];
    ha = [];
    kE.push(kineticE);
    pE.push(potentialE);
    ha.push(hamilton);  
  }
  
    //graph
  if (isDrawing){
    for (let i = 0; i < kE.length; i++) {
      const x = i;
      const y1 = kE[i] * 5;
      const y2 = pE[i] * 5;
      const y3 = ha[i] * 5;
      noStroke()
      fill(150, 150, 0);
      ellipse(x +400, y1 +550, 2, 2);
      fill(0, 150, 150);
      ellipse(x +400, y2 +550, 2, 2);
      fill(150, 0, 150);
      ellipse(x +400, y3 +550, 2, 2);
    }
  
    //diagram
  
  for (let atom = 0; atom < atomCount; atom++) { 
    let x = pos[atom][0] *100
    let y = pos[atom][1] *100
    let z = pos[atom][2] *100
    stroke(0)
    strokeWeight(2);
    fill(myColor[atom])
    ellipse(x, y, 30);
    ellipse(x, z + 340, 30);
    strokeWeight(1);
    fill(255);
    text(`z: ${340-int(Math.floor(z/10) * 10)}`, x, y + 2);
    text(`y: ${340-int(Math.floor(y/10) * 10)}`, x, z + 342);
    // ellipse(z + 340, y, 30);
    }
    
  }else {
    
   for (let atom = 0; atom < atomCount; atom++) { 
    let x = initialPos[atom][0] *100
    let y = initialPos[atom][1] *100
    let z = initialPos[atom][2] *100
    stroke(0)
    strokeWeight(2);
    fill(myColor[atom])
    ellipse(x, y, 30);
    ellipse(x, z + 340, 30);
    strokeWeight(1);
    fill(255);
    text(`z: ${340-int(Math.floor(z/10) * 10)}`, x, y + 2);
    text(`y: ${340-int(Math.floor(y/10) * 10)}`, x, z + 342);
    // ellipse(z + 340, y, 30);
    } 
    
  }
  
    
      // keymap
    let w = 370;
    let l = 280
    rect(w-10, l-35 , 150, 93)
    noStroke();
    fill(150, 0, 0);
    rect(w, l, 50, 50)
    fill(0, 0, 150);
    beginShape();
    vertex(w, l-5);
    vertex(w +50 ,l-5);
    vertex(w +80, l - 25);
    vertex(w +30, l - 25);
    endShape(CLOSE); 
    
    // info
  
  noStroke();
  fill(0)
  textFont();
  
  textAlign(LEFT);
  let s = 40;
  let t = 20;
  let u =  365
  textSize(19);
   text('Sensing Repulsive-Attractive Forces in Molecular Interactions', u, s - t * 1.2);
  textSize(14);
  text('Van der Waals Interaction Modeled by Lennard-Jones Potential', u, s + t * 0);
  text('Incorporating Periodic Boundary Conditions and Minimum Image Convention', u, s + t * 1);
  text(`Frames per Second: ${fps} (dt: 0.002)`, u, s + t * 2);

  text('Temperature: 100K', u, s + t * 4);
  text('Box Size (dimensionless): 6.8', u, s + t * 5);
  text('Cut-off Radius (dimensionless): 2.5', u, s + t * 6);
    // label

  
  text("Kinetic Energy", windowWidth -150, 360);
  text("Potential Energy", windowWidth -150, 380);
  text("Hamiltonian", windowWidth -150, 400);
  textAlign(RIGHT);
  let b = 370
  let c = 360
  text(0, b, 555);
  text(10, b, -50 + 555);
  text(20, b, -100 + 555);
  text(-10, b, 50 + 555);
  text('keymap', w+130, l+49);
  
  
  stroke(0);
  line(c +25, 550, c+30, 550);
  line(c +25, 50 +550, c+30, 50 +550);
  line(c +25, -100 +550, c+30, -100 +550);
  line(c +25, -50 + 550, c+30, -50 + 550);
  line(c +28, -130 + 550, c+28, 100 + 550);
  
  let v1 = createVector( c +28, -130 + 550, c+28, 100 + 550);
  let v2 = createVector( c+28, 100 + 550,  c +28, -150 + 550);
  dimensionArrow(v1, v2);
  dimensionArrow(v2, v1);
  
  stroke(150, 150, 0);
  line(windowWidth-170, 355, windowWidth-160, 355);
  stroke(0, 150, 150);
  line(windowWidth-170, 375, windowWidth-160, 375);
  stroke(150, 0, 150);
  line(windowWidth-170, 395, windowWidth-160, 395);
  
}

function importFile(filename) {
    loadStrings(filename, fileLoaded);
}

function handleFile(file) {
  if (file.type === 'text') {
    loadStrings(file.data, fileLoaded);
  } 
}

function fileLoaded(data) {
  const fileLines = data;
  atomCount = fileLines.length;
  pos = [];
  for (let line = 0; line < atomCount; line++) {
    const coordsList = fileLines[line].split(' ').filter(i => i !== '');
    const coords = coordsList.map(parseFloat);
    pos.push(coords);
  } result = { 
      pos: pos,
      atomCount: atomCount
    
    };
}

function initVel() {
    let vel = [];
    for (let atom = 0; atom < atomCount; atom++) {
        let velX = randomGaussian(0, Math.sqrt(tInit));
        let velY = randomGaussian(0, Math.sqrt(tInit));
        let velZ = randomGaussian(0, Math.sqrt(tInit));
        vel.push([velX, velY, velZ]);
    }
    
    // Adjust velocities to have zero mean
    let meanX = 0, meanY = 0, meanZ = 0;
    for (let atom = 0; atom < atomCount; atom++) {
        meanX += vel[atom][0];
        meanY += vel[atom][1];
        meanZ += vel[atom][2];
    }
    meanX /= atomCount;
    meanY /= atomCount;
    meanZ /= atomCount;
    for (let atom = 0; atom < atomCount; atom++) {
        vel[atom][0] -= meanX;
        vel[atom][1] -= meanY;
        vel[atom][2] -= meanZ;
    }
    return vel;
}

function calculateForces(p) {
    const forces = new Array(atomCount).fill().map(() => [0, 0, 0]);
    let potentialSum = 0;
    let virial = 0;

    for (let atom = 0; atom < atomCount - 1; atom++) {
        for (let other = atom + 1; other < atomCount; other++) {
            const [rx, ry, rz, r] = calculateDist(p, atom, other);
            if (r < rcut) {
                const ucut = ((4 / rcut ** 12) - (4 / rcut ** 6))
                const ducut = -((48 / (rcut ** 13)) - (24 / (rcut ** 7)))
                const denom = Math.pow(r, 6);
                const denomSq = Math.pow(denom, 2);
                
                potentialSum += ((4 / denomSq) - (4 / denom)) - ucut - (r - rcut) * ducut;
                const force = ((48 / (r * denomSq)) - (24 / (r * denom))) + ducut;
                const Fx = force * (rx / r);
                const Fy = force * (ry / r);
                const Fz = force * (rz / r);
                forces[atom][0] -= Fx;
                forces[other][0] += Fx;
                forces[atom][1] -= Fy;
                forces[other][1] += Fy;
                forces[atom][2] -= Fz;
                forces[other][2] += Fz;
                virial += Fx * rx + Fy * ry + Fz * rz;
            }
        }
    }
    return [forces, potentialSum ];
}

function calculateDist(p, atom, other) {
    let rx = p[other][0] - p[atom][0];
    let ry = p[other][1] - p[atom][1];
    let rz = p[other][2] - p[atom][2];

    if (rx <= -boxSize * 0.5) {
        rx += boxSize;
    } else if (rx > boxSize * 0.5) {
        rx -= boxSize;
    }

    if (ry <= -boxSize * 0.5) {
        ry += boxSize;
    } else if (ry > boxSize * 0.5) {
        ry -= boxSize;
    }

    if (rz <= -boxSize * 0.5) {
        rz += boxSize;
    } else if (rz > boxSize * 0.5) {
        rz -= boxSize;
    }

    const r = Math.sqrt(rx ** 2 + ry ** 2 + rz ** 2);
    return [rx, ry, rz, r];
}

function updateVel(v, f) {
    const vNew = new Array(3);
    for (let dir = 0; dir < 3; dir++) {
        vNew[dir] = v[dir] + (dt / 2) * f[dir];
    }
    return vNew;
}

function updatePos(p, v) {
    const pNew = new Array(3);
    for (let dir = 0; dir < 3; dir++) {
        pNew[dir] = p[dir] + dt * v[dir];
        if (pNew[dir] >= boxSize) {
            pNew[dir] -= boxSize;
        } else if (pNew[dir] < 0) {
            pNew[dir] += boxSize;
        }
    }
    return pNew;
}

function calculateMean(velocities) {
    let sum = 0;
    for (let i = 0; i < velocities.length; i++) {
        sum += velocities[i];
    } 
    return sum / velocities.length;
}

function extractVelocities(velocities) {
    let velArray = [];
    for (let i = 0; i < velocities.length; i++) {
        velArray.push(velocities[i][0]);
        velArray.push(velocities[i][1]);
        velArray.push(velocities[i][2]);
    }
    return velArray;
}

function calculateStandardDeviation(array, mean) {
    const variances = array.map(val => Math.pow(val - mean, 2));
    const varianceSum = variances.reduce((acc, val) => acc + val, 0);
    return Math.sqrt(varianceSum / array.length);
}

function dimensionArrow(base, vec) {
  push();
  fill(0);
  vec1 = p5.Vector.sub(vec, base);
  translate(base.x, base.y);
  rotate(vec1.heading());
  let arrowSize = 6;
  translate(vec1.mag() - arrowSize, 0); //move vector direction?
  triangle(0, arrowSize / 3, 0, -arrowSize / 3, arrowSize, 0);
  pop();
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
}

function startDrawing() {
  isDrawing = true;
}

function stopDrawing() {
  isDrawing = false;
}