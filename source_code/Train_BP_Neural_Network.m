function [weight1,weight2,bias1,bias2,loss]=Train_BP_Neural_Network(input,weight1,weight2,bias1,bias2,label,lr)
    sample_num=size(input,1);
    input_neuro_num=size(input,2);
    neuro_num_1=size(weight1,2);
    neuro_num_2=size(weight2,2);
    for i=1:sample_num
        Z1=input(i,:)*weight1+bias1;
        A1=sigmoid(Z1);
        Z2=A1*weight2+bias2;
        A2=sigmoid(Z2);
        [loss,d_A2]=CrossEntropyLoss(label(i,:),A2);
        d_Z2=d_A2'.*d_sigmoid(Z2);
        d_weight2=zeros(neuro_num_1,neuro_num_2);
        for j=1:neuro_num_2
            d_weight2(:,j)=(d_Z2(j)*A1)';
            d_bias2=d_Z2;
        end
        d_A1=d_Z2*(weight2');
        d_Z1=d_A1.*d_sigmoid(Z1);
        d_weight1=zeros(input_neuro_num,neuro_num_1);
        for j=1:neuro_num_1
            d_weight1(:,j)=(d_Z1(j)*input(i,:))';
            d_bias1=d_Z1;
        end
        weight1=weight1-lr*d_weight1;
        bias1=bias1-lr*d_bias1;
        weight2=weight2-lr*d_weight2;
        bias2=bias2-lr*d_bias2;
    end
end

function output=d_sigmoid(input)
    output=sigmoid(input).*(1-sigmoid(input));
end

function output=sigmoid(input)
    output=1./(1.+exp(-1.*input));
end

