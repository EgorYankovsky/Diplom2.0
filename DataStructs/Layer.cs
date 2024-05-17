namespace DataStructs;


public record Layer(double Z0, double Z1, double Mu, double Sigma);

public record Anomaly(double X0, double X1, double Y0, double Y1, double Z0, double Z1, double Mu, double Sigma);