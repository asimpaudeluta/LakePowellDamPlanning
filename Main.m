showPower();
function showPower()
    synthData = readmatrix('Project/PowerOutputAR.csv'); %AutoRegressive
    %synthData = readmatrix('Project/PowerOutput.csv');
    Jan = synthData(3,2:601);
    Feb = synthData(4,2:601);
    Mar = synthData(5,2:601);
    Apr = synthData(6,2:601);
    May = synthData(7,2:601);
    Jun = synthData(8,2:601);
    Jul = synthData(9,2:601);
    Aug = synthData(10,2:601);
    Sep = synthData(11,2:601);
    Oct = synthData(12,2:601);
    Nov = synthData(13,2:601);
    Dec = synthData(14,2:601);
    figure();
    ecdf(Jan);
    hold on
    ecdf(Feb);
    hold on
    ecdf(Mar);
    hold on
    ecdf(Apr);
    hold on
    ecdf(May);
    hold on
    ecdf(Jun);
    hold on
    ecdf(Jul);
    hold on
    ecdf(Aug);
    hold on
    ecdf(Sep);
    hold on
    ecdf(Oct);
    hold on
    ecdf(Nov);
    hold on
    ecdf(Dec);
    legend('January','February','March','April','May','June','July','August','September','October','November','December','Location','southeast');
end
function answer = PowerProblem()
    synthData = readmatrix('Project/SynthAR.csv');
    %OptimumRelease checks
    Or = [1417700,965560,783180,784850,1227900,100000,100000,100000,700000,350000,350000,2355000]; %Change this for every situation based on answer from question 1
    E = [321.137969765971,376.88426875,579.788981024669,949.008423529413,1121.58932890576,1911.87097124183,2272.56571853257,2228.06223845667,2078.71896732026,1381.85365654649,1361.42203071895,1039.66365338393];
    answer =  zeros(12,600);
    for j = 2:30:601
        %Initialize each set
        index = 1;
        P(index) = 0;
        St(index) = 0.8*25160000;%Starts from 0, just constructed
        K = 25160000; %Capacity
        for i= j:j+29
            %Initialize each year
            for k= 1:12 %Each month
                %disp(index+":"+i+","+(k+2)+","+synthData(k+2,i));
                Q = getAcrFt(synthData(k+2,i)); %Inflow every year,month
                checkVar = St(index)+Q-E(k)-Or(k); %Level in reservoir
                Release = 0;
                %Release Decision
                if checkVar>K
                    Release = St(index)+Q-E(k)-K;
                elseif checkVar<=K && checkVar>=0
                    Release = Or(k);
                else
                    Release = St(index)+Q-E(k);
                end
                St(index+1) = St(index)+Q-E(k)-Release;
                P(index) = 0.65*1000*9.81*toFlowSI(Release)*hFromUSGS(St(index));
                answer(k,i-1) = P(index);
                index = index + 1;
            end
        end
    end
    %disp(answer)
end
function checkSyntheticData()
    %synthData = readmatrix('Project/syntheticData.csv');
    synthData = readmatrix('Project/SynthAR.csv'); %Auto Regression
    %OptimumRelease checks
    Or = [1417700,965560,783180,784850,1227900,100000,100000,100000,700000,350000,350000,2355000]; %Change this for every situation based on answer from question 1
    E = [321.137969765971,376.88426875,579.788981024669,949.008423529413,1121.58932890576,1911.87097124183,2272.56571853257,2228.06223845667,2078.71896732026,1381.85365654649,1361.42203071895,1039.66365338393];
    NoOfLessThanOrForYear = zeros(1,30);
    for j = 2:30:601
        %Initialize each set
        index = 1;
        St(index) = 0;%Starts from 0, just constructed
        K = 25160000; %Capacity
        for i= j:j+29
            %Initialize each year
            NoOfLessThanOr = 0;
            for k= 1:12 %Each month
                %disp(index+":"+i+","+(k+2)+","+synthData(k+2,i));
                Q = getAcrFt(synthData(k+2,i)); %Inflow every year,month
                checkVar = St(index)+Q-E(k)-Or(k); %Level in reservoir
                Release = 0;
                %Release Decision
                if checkVar>K
                    Release = St(index)+Q-E(k)-K;
                elseif checkVar<=K && checkVar>=0
                    Release = Or(k);
                else
                    Release = St(index)+Q-E(k);
                end
                St(index+1) = St(index)+Q-E(k)-Release;
                if Release<0.8*Or(k) % 80% of required release, put 1 when not required
                    NoOfLessThanOr = NoOfLessThanOr + 1;
                    % disp("Found one Suboptimal");
                end
                index = index + 1;
                % disp(Release-Or(k));
            end
            %disp(i-j+1) %0 to 30
            NoOfLessThanOrForYear(i-j+1) = NoOfLessThanOrForYear(i-j+1)+NoOfLessThanOr;
            % disp("_end of a year");
        end
        % disp(1+((j-2)/30))
        % disp("_end of a set");
    end
    %disp(NoOfLessThanOrForYear);
    answer = "";
    for i = 1:30
        answer = answer + NoOfLessThanOrForYear(i)+",";
    end
    disp(answer);
end
function showSynthData()
    synthData = readmatrix('Project/SynthAR.csv'); %autoregression
    %synthData = readmatrix('Project/syntheticData.csv');
    Jan = synthData(3,2:601);
    Feb = synthData(4,2:601);
    Mar = synthData(5,2:601);
    Apr = synthData(6,2:601);
    May = synthData(7,2:601);
    Jun = synthData(8,2:601);
    Jul = synthData(9,2:601);
    Aug = synthData(10,2:601);
    Sep = synthData(11,2:601);
    Oct = synthData(12,2:601);
    Nov = synthData(13,2:601);
    Dec = synthData(14,2:601);
    figure();
    ecdf(Jan);
    hold on
    ecdf(Feb);
    hold on
    ecdf(Mar);
    hold on
    ecdf(Apr);
    hold on
    ecdf(May);
    hold on
    ecdf(Jun);
    hold on
    ecdf(Jul);
    hold on
    ecdf(Aug);
    hold on
    ecdf(Sep);
    hold on
    ecdf(Oct);
    hold on
    ecdf(Nov);
    hold on
    ecdf(Dec);
    legend('January','February','March','April','May','June','July','August','September','October','November','December','Location','southeast');
end
function genearateSyntheticData()
    %reading csv file flows and generating best distributions and simple synthetic
    flowsCsv = readmatrix('Project/DailyInflowSortedbyMonth.csv');
    FittedRand = "";
    for i = 2:13
        [dist,p,e40,e50,e60,m1,m2,m3,fR] = getOutputA(i,flowsCsv);
        %disp(dist+","+p+","+e40+","+e50+","+e60+","+m1+","+m2+","+m3);
         for j = 1:600
             FittedRand = FittedRand+fR(j)+",";
         end
         FittedRand = FittedRand+char(13);
    end
         disp(FittedRand)
    function [MyDistribution,pmax,e40,e50,e60,m1,m2,m3,fR] = getOutputA(i,inputFile)
        [pN,e40N,e50N,e60N,m1N,m2N,m3N,fRN] = distributionFunction(inputFile(:,i),"Normal");
        [pLn,e40Ln,e50Ln,e60Ln,m1Ln,m2Ln,m3Ln,fRLn] = distributionFunction(inputFile(:,i),"LogNormal");
        [pGev,e40Gev,e50Gev,e60Gev,m1Gev,m2Gev,m3Gev,fRGev] = distributionFunction(inputFile(:,i),"gev");
        [pGam,e40Gam,e50Gam,e60Gam,m1Gam,m2Gam,m3Gam,fRGam] = distributionFunction(inputFile(:,i),"Gamma");
        MyDistribution = "Normal";
        pmax = pN;
        e40 = e40N;
        e50 = e50N;
        e60 = e60N;
        m1 = m1N;
        m2 = m2N;
        m3 = m3N;
        fR = fRN;
        if pLn>pmax
            MyDistribution = "LogNormal";
            pmax = pLn;
            e40 = e40Ln;
            e50 = e50Ln;
            e60 = e60Ln;
            m1 = m1Ln;
            m2 = m2Ln;
            m3 = m3Ln;
            fR = fRLn;
        end
        if pGev>pmax
            MyDistribution = "gev";
            pmax = pGev;
            e40 = e40Gev;
            e50 = e50Gev;
            e60 = e60Gev;
            m1 = m1Gev;
            m2 = m2Gev;
            m3 = m3Gev;
            fR = fRGev;
        end
        if pGam>pmax
            MyDistribution = "Gamma";
            pmax = pGam;
            e40 = e40Gam;
            e50 = e50Gam;
            e60 = e60Gam;
            m1 = m1Gam;
            m2 = m2Gam;
            m3 = m3Gam;
            fR = fRGam;
        end
    function [p,e40,e50,e60,m1,m2,m3,fittedRandom] = distributionFunction(inputData,distribution)
            pd = fitdist(inputData,distribution);
            xvals = 0:max(inputData)/100:max(inputData);
            y = cdf(pd,xvals);
            % figure();
            % plot(xvals,y);
            % hold on;
            % ecdf(inputData);
            % legend(distribution,'Empirical','Location','northwest');
            % figure();
            % qqplot(inputData,pd);
            [~,p,~,~] = kstest(inputData,'CDF',pd);
            e40 = icdf(pd,0.4);
            e50 = icdf(pd,0.5);
            e60 = icdf(pd,0.6);
            m1 = mean(pd);
            m2 = var(pd);
            m3 = skewness(inputData);
            %synthetic data simple
            matrixRandomUniform = rand(600);
            uniformDist = matrixRandomUniform(:,1);
            fittedRandom = icdf(pd,uniformDist);
        end
    end
end
function storageOptimProblem()
    %one c
    Ft = 0.7 * 25160000;
    Rect = 0.9 * 25160000;
    Rt1 = 350000;
    Rt2 = 2*Rt1;
    outputFile = readtable("Project/Output.csv");
    flows1c = outputFile{:,5};
    evap1c = outputFile{1,23:34};
    release1c = outputFile{1,35:46}; %unused    
    st(1)=0.8*25160000;
    for i = 0:35
        %decision variables
        release = optimvar('release','LowerBound',100000);
        storage = optimvar('storage','LowerBound',0,'UpperBound',25160000);
        %objective function
        Wf = 0; Wr1 = 0;  Wr2 = 0; Wrec = 0;
        switch mod(i,12)
            case 0 %January
            Wf = 1;
            Wr1 = 1;
            case 1 %February
            Wf = 1;
            Wr1 = 1;
            case 2 %March
            Wf = 1;
            Wr1 = 1;
            case 3 %April
            Wf = 1;
            Wr1 = 1;
            case 4 %May
            Wf = 1;
            Wr1 = 1;
            case 5 %June
            Wr2 = 1;
            Wrec = 1;
            case 6 %July
            Wr2 = 1;
            Wrec = 1;
            case 7 %August
            Wr2 = 1;
            Wrec = 1;
            case 8 %September
            Wr2 = 1;
            case 9 %October
            Wr1 = 1;
            case 10 %November
            Wr1 = 1;
            case 11 %December
            Wr1 = 1;
            Wf = 1;
        end
        
        inflow = getAcrFt(flows1c(mod(i,12)+1));
        loss = evap1c(mod(i,12)+1);
        prob = optimproblem('Objective',Wf*(myMax(st(i+1)-Ft,0))^2+Wf*(myMax(storage-Ft,0))^2+Wrec*(myMin(storage-Rect,0))^2+Wrec*(myMin(st(i+1)-Rect,0))^2+Wr1*(myMin(release-Rt1,0))^2+Wr2*(myMin(release-Rt2,0))^2,'ObjectiveSense','min');
        %Constraints
        prob.Constraints.cons1 = storage == st(i+1) + inflow - loss - release;
        %Solving
        problem = prob2struct(prob);
        [sol,~,~,~] = quadprog(problem);
        %change storage for next iteration
        ReleaseValues(i+2) = sol(1);
        st(i+2) = sol(2);
    end
    
    function max = myMax(a,b)
        try
           if double(a)<=b
                max = b;
            else
                max = a;
           end
           disp("Working")
        catch
           max = a;
        end   
    end
    function min = myMin(a,b)
        try
            if double(a)>=b
                min = b;
            else
                min = a;
            end
        catch
            min = a;
        end
    end
end
function acrFt = getAcrFt(cfs)
    acrFt = cfs*(30*86400)/43560;
end
function SIFlow = toFlowSI(acrFt)
    SIFlow = acrFt * 0.00047587962;
end
function generateSyntheticDataAutoRegression()
    %synthetic data lag 1 autocorrelation
    inputFile = readmatrix('Project/inputAR.csv');
    x = inputFile(2,16); %mu
    s = inputFile(2,17); %sigma
    p = inputFile(2,18); %rho
    normalPd = makedist("Normal","mu",0,"sigma",1);
    pd_InflowNormal = fitdist(inputFile(6:56,2),"Normal");
    for index = 1:30:600
        %Every set
        Q(index) = icdf(pd_InflowNormal,rand(1));
        for each = index:index+29
            %Every year
            N0_1(each) = icdf(normalPd,rand(1));
            Q1 = Q(each);
            Q2 = x+p*(Q1-x)+N0_1(each)*(((s^2)*(1-(p^2)))^0.5);
            Q(each + 1) = Q2;
        end
    end
    for i = 1:600
        for j = 1:12
            mon(j) = sepFrac(j)*Q(i);
        end
        disp(Q(i)+","+mon(1)+","+mon(2)+","+mon(3)+","+mon(4)+","+mon(5)+","+mon(6)+","+mon(7)+","+mon(8)+","+mon(9)+","+mon(10)+","+mon(11)+","+mon(12));
    end
    function frac = sepFrac(month)
        frac = inputFile(3,month+2)*365/inputFile(60,month+2);
    end
plot(Q)
end
function h = hFromUSGS(storageAcreFt)
    h = 36.592*log(storageAcreFt)-476.92;
end
