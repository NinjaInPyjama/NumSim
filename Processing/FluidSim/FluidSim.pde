import java.io.*;

final float width = 1000;
final float height = 1000;

Geometry geom;
Parameter param;
Compute comp;

int gridID;
Visu visu;

void setup() {
  size(1000,1000);
  
  geom = new Geometry("default_geom");
  param = new Parameter("default_param");
  comp = new Compute(geom, param);
  
  gridID = 4;
  
  visu = new Visu(width, height, geom);
  visu.setResolution(new int[]{200, 200});
  visu.setGrid(comp.getVelocity());
  colorMode(HSB, 360);
  
  MultiGrid mg = new MultiGrid(geom);
  mg.printCellSize();
  mg = mg.restrict(1);
  mg.printCellSize();
  mg = mg.restrict(2);
  mg.printCellSize();
  mg = mg.restrict(4);
  mg.printCellSize();
  mg = mg.restrict(8);
  mg.printCellSize();
  mg = mg.interpolate(16);
  mg.printCellSize();
  mg = mg.interpolate(8);
  mg.printCellSize();
  mg = mg.interpolate(4);
  mg.printCellSize();
  mg = mg.interpolate(2);
  mg.printCellSize();
}

void draw() {
  background(0,0,0);
  
  if(comp.getTime() < param.getTend()) {
    comp.timeStep(false);
    
    switch(gridID) {
      case 1: 
        visu.setGrid(comp.getU().copy());
        break;
      case 2: 
        visu.setGrid(comp.getV().copy());
        break;
      case 3: 
        visu.setGrid(comp.getP().copy());
        break;
      case 4: 
        visu.setGrid(comp.getVelocity().copy());
        break;
      case 5: 
        visu.setGrid(comp.getVelocity().copy());
        visu.setPathLine(comp.getPathLines());
        visu.showPathLines(true);
        break;
      case 6: 
        visu.setGrid(comp.getVelocity().copy());
        visu.setStreakLine(comp.getStreakLines());
        visu.showStreakLines(true);
        break;
      case 7: 
        visu.setGrid(comp.getVorticity().copy());
        break;
      case 8: 
        visu.setGrid(comp.getStream());
        break;
      case 9: 
        visu.setGrid(comp.getStream());
        break;
      default: 
        visu.setGrid(comp.getVelocity().copy());
        break;
    }
    
    visu.draw();
    
    fill(0);
    textAlign(LEFT, TOP);
    textSize(24);
    //text("Simulated Time: " + comp.getTime() + "s", 10, 10);
    //text("Max Value: " + visu._maxValue, 10, 40);
  }
}

void keyPressed() {
  visu.removeNiveauLines();
  visu.showPathLines(false);
  visu.showStreakLines(false);
  
  switch(key) {
    case '1':
      gridID = 1;
      break;
    case '2':
      gridID = 2;
      break;
    case '3':
      gridID = 3;
      break;
    case '4':
      gridID = 4;
      break;
    case '5':
      gridID = 5;
      break;
    case '6':
      gridID = 6;
      break;
    case '7':
      gridID = 5;
      break;
    case '8':
      gridID = 6;
      break;
    case '9':
      gridID = 6;
      visu.showNiveauLines(50);
      break;
  }
}