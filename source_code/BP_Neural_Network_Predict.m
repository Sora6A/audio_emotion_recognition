function output=BP_Neural_Network_Predict(weight1,weight2,bias1,bias2,input)
    sample_num=size(input,1);
    neuro_num_2=size(weight2,2);
    output=zeros(sample_num,neuro_num_2);
    for i=1:sample_num
        Z1=input(i,:)*weight1+bias1;
        A1=sigmoid(Z1);
        Z2=A1*weight2+bias2;
        A2=sigmoid(Z2);
        output(i,:)=A2;
    end
end

function output=sigmoid(input)
    output=1./(1.+exp(-1.*input));
end
