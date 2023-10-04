package odesolver.exception;

public class InterpolationException extends RuntimeException{
    private String errorMsg;

    public InterpolationException(int itask, double tout){
        super();
        errorMsg = "lsoda: trouble from intdy, itask = "+itask+
                ", tout = "+tout;
    }

    @Override
    public String getMessage() {
        return errorMsg;
    }
}
