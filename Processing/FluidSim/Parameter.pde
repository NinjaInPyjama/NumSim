class Parameter {
  private float _re;
  private float _alpha;
  private float _dt;
  private float _tend;
  private int _itermax;
  private float _eps;
  
  private final String SAVE_FILE_LOCATION;
  
  public Parameter(String filename) {
    SAVE_FILE_LOCATION = "data/" + filename + ".json";
    load();
  }
  
  private void load() {    
    JSONObject json = loadJSONObject(SAVE_FILE_LOCATION);
    _re = json.getFloat("re");
    _alpha = json.getFloat("alpha");
    _dt = json.getFloat("dt");
    _tend = json.getFloat("tend");
    _itermax = json.getInt("itermax");
    _eps = json.getFloat("eps");
  }
  
  public float getRe() {
    return _re;
  }
  
  public void setRe(float re) {
    _re = re;
  }
  
  public float getAlpha() {
    return _alpha;
  }
  
  public float getDt() {
    return _dt;
  }
  
  public float getTend() {
    return _tend;
  }
  
  public int getItermax() {
    return _itermax;
  }
  
  public float getEps() {
    return _eps;
  }
}