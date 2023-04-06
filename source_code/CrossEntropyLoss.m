function [loss,dloss]=CrossEntropyLoss(label,predict_output)
    loss=0;
    for i=1:length(label)
        loss=loss-label(i)*log(predict_output(i))-(1-label(i))*log(1-predict_output(i));
    end
    dloss=zeros(length(label),1);
    for i=1:length(label)
        dloss(i)=-label(i)/predict_output(i)+(1-label(i))/(1-predict_output(i));
    end
end