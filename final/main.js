var canvas1={};

canvas1.setup = function() {
    
    this.setGlobals();
    
    this.cthulhu = true;
    this.uni =     true;
    
    if(this.uni)
        this.makeUni();
    if(this.cthulhu)
        this.makeCthulhu();
    
    //using these for eyes
    this.makeTestObjects();
    
    
    // leaf.setSampler2D("uSampler", "img/maple_diffuse.jpg");
    // leaf.setSampler2D("uAlphaSampler", "img/maple_alpha.jpg");
}

canvas1.update = function(){
    
    ////***Setup***////
    
    this.noise = (this.n.noise(time/2,2)*.01);
    
    var mouseX = this.mouseX/512-1,
    mouseY = -this.mouseY/512+1;
    
    this.mousex = mouseX;
    this.mousey = mouseY;
    
    var o = this.things;
    
    if(mouseX - this.pMouse[0]>0 || mouseX - this.pMouse[0]<0){
        this.speed += (mouseX-this.pMouse[0])*PI;
    }
    
    //add current location to average
    this.a.add((mouseX-this.pMouse[0])*PI);
    
    ///**__**__**__**/Translate Objects/**__**__**__**////
    if(this.uni)
        this.moveUni(o,this.speed,mouseX,mouseY);
    if(this.cthulhu){
        this.moveCthulhu(mouseX,mouseY);
    }
    
    //move environment
    var m = 15.3;
    o.torusParent.move(o.mountains,makeMatrix({z:-9,rx:PI/2}));
    o.mountains.tMatrix = (makeMatrix({sx:m,sy:m,sz:m,z:.343,rz:-.1*this.rotSpeed}));
    //o.mountains.tMatrix = (makeMatrix({sx:m,sy:m,sz:m,z:.343,rz:time*.2}));
    this.things.mountains.setUniform('timer',time*.05);
    this.things.mountains.setUniform('mouser',-this.mousex);
    
    //this.things.drop.move(this.things.drop,makeMatrix({rx:mouseX,ry:mouseY})); //move drop
    //this.particles(10000); not working that well and slow
    
    /////**__**__**__**////Management/**__**__**__**////
    
    this.pMouse.push(mouseX+this.noise);
    this.pMouse.push(mouseY);
    
    while(this.pMouse.length>2){
        this.pMouse.shift();
    }
    if(this.pMouse[0]==mouseX){
        this.pMouse[0]=mouseX;
        this.pMouse[1]=mouseY;
    }
}

canvas1.makeSpringyTentacles = function(num,joints){
    for(var i = 1 ; i < num; i++){
        //pTents has no geometry - just placeholders
        this.addTent(this.pTents,joints,'none',i);
        if(i<7)
            this.addTent(this.spTents,joints,'sphere',i,6);
        else
            this.addTent(this.spTents,joints,'drop',i,12);
        
    }
    for(var i = 0 ; i < this.spTents.length ; i++){
        for(var j = 0 ; j < this.spTents[i].length ; j++){
            var s = this.pTents[i].length
            this.spTents[i][j].spring = new Spring({P:[1,0,0],mass:map_range(j,1,s,.2,3),damp:Math.pow(map_range(j,1,s,3,.5),.45)});
        }
    }
}

canvas1.addTent = function(array,num,typer,id,div,parms){
    
    var type = typer || 'cube';
    var d = div || 12;
    var tent = [];
    
    for(var i = 0 ; i < num ; i++){
        
        if(type=='cube')
            this.addObject(createCube(), 'fs_phong');
        else if(type=='none')
            this.addObject([], 'fs_phong');
        else if(type=='drop')
            this.addObject(createParametric(1/(d*2),1/d,drop,{r1:map_range(i,0,num,.01,.2)}), 'squishy2');
        else
            this.addObject(createParametric(1/(d*2),1/d,sph,{r1:map_range(i,0,num,.05,.01)}), 'squishy');
        
        tent.push(this.objects[this.objects.length-1]);
        this.things["tent"+id+i]=this.objects[this.objects.length-1];
        
        //sadly the uniforms are hardcoded - it would be great if I could pass them in an argument
        if(type=='drop'){
            this.objects[this.objects.length-1].setUniform('p', [.004,0.008,0.0041,0,0,0,0,0,0.001,i+2]);
            for(var q = 1 ; q < this.objects[this.objects.length-1].vertices.length ; q+=8){
                this.objects[this.objects.length-1].vertices[q] *= map_range(i,0,num,8,1);
            }
            this.objects[this.objects.length-1].vertexBuffer = createVertexBuffer(this.canvas.gl, this.objects[this.objects.length-1].vertices);
            
        }
        else
            this.objects[this.objects.length-1].setUniform('p', [-.05+(i/num*.15),-.05+(i/num*.12), .002,.01,.02,.03,.13,.2,.03,i+2]);
        
        this.objects[this.objects.length-1].setUniform('lDir',[1,.3,-.5]);
        this.objects[this.objects.length-1].setUniform('off',i*1.0);
        if(i>0){
            tent[i-1].add(tent[i]);
        }
    }
    
    array.push(tent);
}

canvas1.makeTestObjects = function(){
    
    this.prti = {
    px:[],
    py:[],
    pz:[],
    vx:[],
    vy:[],
    vz:[]
    };
    //testDrop
    /*
     this.makeThing(
     {rx:Math.PI/6},
     "drop",
     [.5,.1, .5,0.5,0.5,0.5, 0.1,0.1,0.1,1],
     createParametric(1/30,1/30,drop,{r1:.41,r2:.5,mp:1,ap:.5}),
     'squishy'
     );
     this.things.drop.matrix = this.things.drop.tMatrix;
     
     this.addObject(createParametric(1/(d*2),1/d,drop,{r1:map_range(i,0,num,.01,.2)}), 'squishy');
     this.objects[this.objects.length-1].setUniform('p', [.01,.008, .002,.01,.02,.03,.13,.2,.03,1]);
     
     
     
     this.makeThing(
     {rx:Math.PI/6},
     "drop",
     [.01,.008, .002,.01,.05,.03,.13,.2,.13,1],
     createParametric(1/30,1/30,sph,{r1:.41,r2:.5,mp:1,ap:.5}),
     'squishy'
     );
     this.things.drop.matrix = this.things.drop.tMatrix;
     this.things.drop.setUniform('lDir',[1,0,-.3]);
     */
    
    this.n3d=new SimplexNoise();
    
    //environment
    this.makeThing(
                   {rx:PI/2},
                   "mountains",
                   [.05,.01, .05,0.05,0.05,0.05, 0.1,0.1,0.1,1],
                   createParametric(1/190,1/190,sph,{r1:1,r2:.5,mp:1,ap:.5}),
                   'texture' //should be 'texture'
                   );
    /*
     this.makeThing(
     {rx:PI/2},
     "sky",
     [.05,.01, .05,0.05,0.05,0.05, 0.1,0.1,0.1,1],
     [],
     'squishy'
     );
     */
    var min = 0;
    var max = 0;
    f.overflow=1000;
    this.things.mountains.verts2 = [];
    this.things.mountains.verts3 = [];
    
    for(var v = 2 ; v < this.things.mountains.vertices.length ; v+= 8){
        var p = this.things.mountains.vertices;
        
        var pusher;
        var nscale = Math.max(0,map_range(p[v],-1,-.7,.02,-.020));
        var nscale2 = Math.pow(Math.max(0,map_range(p[v],-1,-.99,.31,0)),2);
        
        if(p[v]<-.7)
            pusher=true;
        
        if(p[v]<min)
            min=p[v];
        if(p[v]>max)
            max=p[v];
        
        var size = .4
        
        if(this.things.mountains.vertices[v]<-size)
            this.things.mountains.vertices[v] = -size+nscale2+(nscale*((this.n3d.noise3d(p[v-2]*6,p[v-1]*6,p[v]*8))+(this.n3d.noise3d(p[v-2]*16,p[v-1]*16,p[v]*16))));
        if(this.things.mountains.vertices[v]>size)
            this.things.mountains.vertices[v] = size;
        
        
        if(pusher){
            for(var q = -2 ; q < 6 ; q++){
                this.things.mountains.verts2.push(p[v+q]);
                pusher=false;
            }
        }
        else
            this.things.mountains.verts3.push(p[v+q]);
    }
    
    this.things.mountains.vertexBuffer = createVertexBuffer(this.canvas.gl, this.things.mountains.verts2);
    this.things.mountains.textures = [];
    this.things.mountains.textures[0] = "images/soil.jpg";
    this.things.mountains.textures[1] = "images/cracks.jpg";
    this.things.mountains.textures[2] = "images/splat2.jpg";
    this.things.mountains.fov =2;
    
    //this.things.sky.vertexBuffer = createVertexBuffer(this.canvas.gl, this.things.mountains.verts3);
    
    this.makeThing(
                   {},
                   "torusParent",
                   [.5,.1, .5,0.5,0.5,0.5, 0.1,0.1,0.1,1],
                   [],
                   'squishy'
                   );
    this.things.torusParent.add(this.things.mountains);
    // this.things.torusParent.add(this.things.sky);
    
    this.makeThing(
                   {rx:Math.PI,z:-4,sz:.4},
                   "sky",
                   [1,.0, .0,0,0,0, 0,0,0,1],
                   createParametric(1/30,1/18,sph,{r1:5,r2:.2}),
                   'flatTex' //should be 'flatTex'
                   );
    
    this.things.sky.textures = [];
    this.things.sky.fov =1;
    this.things.sky.textures[0] = "images/sky.jpg";
    this.things.sky.matrix=this.things.sky.tMatrix;
    this.things.mountains.add(this.things.sky);
    
    var particlePusher = [];
    for(var i = 0 ; i <100 ; i++){
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-7]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-6]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-5]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-4]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-3]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-2]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length-1]);
        particlePusher.push(this.things.mountains.verts2[this.things.mountains.verts2.length]);
    }
    
    this.makeThing(
                   {},
                   "particles",
                   [1,.0, .0,0,0,0, 0,0,0,1],
                   particlePusher,
                   'fs_phong'
                   );
    
    
    this.things.particles.drawType="POINTS";
    
    //eyeballs
    
    this.things.mountains.setUniform('lDir',[0,1,0]);
    //this.things.mountains.setUniform('mouser',this.mousex);
    
    this.makeThing(
                   {},
                   "testBall1",
                   [.0,.0, .0,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,sph,{r1:.021,r2:.2}),
                   'fs_phong'
                   );
    this.makeThing(
                   {},
                   "testBall2",
                   [1,0, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,sph,{r1:.021,r2:.2}),
                   'eye'
                   );
    this.things.testBall2.textures = [];
    this.things.testBall2.textures[0] = "images/eye.jpg";
    this.things.testBall2.textures[1] = "images/beach.jpg";
}

canvas1.setGlobals = function(){
    this.cthulhu = true;
    this.uni = true;
    //springs for body
    this.sp  = new Spring({P:[1,0,0],mass:3,damp:.9});
    this.sp2 = new Spring({P:[1,0,0],mass:3,damp:.9});
    //things is for hashing
    this.things = [];
    //previous mouse position
    this.pMouse = [0,0];
    this.speed = 0;
    //this.loff = this.roff = this.luxoff = this.luyoff = this.ruxoff = this.ruyoff = 0;
    this.n = new SimplexNoise();
    this.a = new Average(10);
    this.tents = [];
    //parent tentacles
    this.pTents = [];
    //tentacles with springs
    this.spTents = [];
    
    this.addSpeed = 3;
    
    this.rotSpeed = 0;
    
    this.stick = new Array(4);
    this.noise;
    this.mnoise=function(t,s,b){
        return map_range(this.n.noise(t,t),0,1,s,b);
        
    }
}

canvas1.makeCthulhu = function(){
    
    //make 3 tentacles and give them their uniforms
    
    for (var p = 0 ; p < 3; p++){
        this.stick[p] = [0,0,0];
        var tent = [];
        var num = 30;
        if(p>1) num = 20;
        if(p>2) num = 20;
        
        for (var i = 0 ; i < num ; i++){
            this.addObject(createParametric(1/(48*2),1/48,sph), 'squishy');
            this.objects[this.objects.length-1].name = "ball" + p + i;
            this.things[this.objects[this.objects.length-1].name]=this.objects[this.objects.length-1];
            this.objects[this.objects.length-1].setUniform('p', [.01,.008, .002,.01,.02,.03, .13,.2,.03,i]);
            this.objects[this.objects.length-1].setUniform('lDir',[1,.3,-.5]);
            this.objects[this.objects.length-1].setUniform('off',i*1.0);
            
            tent.push(this.objects[this.objects.length-1]);
            
        }
        this.tents.push(tent);
    }
    
    this.makeSpringyTentacles(9,30);
}

canvas1.moveUni = function(o,speed,mouseX,mouseY){
    
    //wheel is the overall parent
    o.wheel.tMatrix=makeMatrix({sx:.25,sy:.25,sz:.25})
    
    var params={rx:PI/2,ry:time};
    var landSpeed = time*this.addSpeed;
    
    //outerTire controls the spinning of the wheel
    o.outerTire.tMatrix=(makeMatrix({rz:speed+landSpeed}));
    o.pedalJoint1.tMatrix=(makeMatrix({x:-.62,z:.92,rz:-speed-landSpeed}));
    o.pedalJoint2.tMatrix=(makeMatrix({x:.62,z:-.92,rz:-speed-landSpeed}));
    
    var a = this.noise;
    
    obj.move(o.wheel,makeMatrix({
                                x:mouseX+a,
                                y:-.6,
                                ry:(sin(landSpeed)*.2)+sin(this.a.avg()*6+mouseX*this.a.avg())
                                }));
    
    this.rotSpeed += a/1.5;
    this.rotSpeed -= (this.addSpeed/3)*(.015);
    
    
    //the wangles are the shaft that go up to the seat
    o.wangleParent.tMatrix=(jMatMult(o.wheel.tMatrix,makeMatrix({rz:this.a.avg()*-3,sy:4,sx:3.5,sz:3.5})));
    o.wangleton.move(o.wangleton,o.wangleton.tMatrix);
}

canvas1.makeThing = function(parms,name,uniform,obj,material){
    this.addObject(obj, material);
    var obj = this.objects[this.objects.length-1];
    obj.name = name;
    obj.tMatrix=makeMatrix(parms);
    obj.setUniform('p', uniform);
    this.things[name] = obj;
}

canvas1.moveCthulhu = function(mouseX,mouseY){
    
    var o = this.things;
    var D = [mouseY,mouseX,0];
    var off = -.1;
    var offy = .8;
    
    var lstick = 0, rstick = 0;
    
    //animate legs
    C = [o.seatPoint.matrix[12],o.seatPoint.matrix[13]-.05,o.seatPoint.matrix[14]];
    R = [ o.pedalJoint1.matrix[12],
         o.pedalJoint1.matrix[13]+.01,
         o.pedalJoint1.matrix[14]];
    var Q = [mouseX+.5,-.2,.2];
    
    this.moveTentacle(0,C,R,Q,".2+(j+1)/25");
    
    C = [o.seatPoint.matrix[12],o.seatPoint.matrix[13]-.05,o.seatPoint.matrix[14]];
    R = [ o.pedalJoint2.matrix[12],
         o.pedalJoint2.matrix[13]+.01,
         o.pedalJoint2.matrix[14]];
    Q = [mouseX+.5,-.2,-.2];
    
    this.moveTentacle(1,C,R,Q,".2+(j+1)/25");
    
    //animate Head
    this.sp.F[0]=o.seatPoint.matrix[12]+.3||0;
    this.sp.F[1]=o.seatPoint.matrix[13]+.2||0;
    this.sp.update(.2);
    
    this.sp2.F[0]=this.sp.P[0]-.5||0;
    this.sp2.F[1]=this.sp.P[1]+.2||0;
    this.sp2.update(.1);
    
    o.testBall1.matrix=makeMatrix({x:o.ball219.matrix[12]-.0,y:o.ball219.matrix[13]+.04,z:.3});
    o.testBall2.matrix=makeMatrix({x:o.ball219.matrix[12]-.01,y:o.ball219.matrix[13]+.03,z:.3});
    
    R = [o.seatPoint.matrix[12],o.seatPoint.matrix[13],o.seatPoint.matrix[14]];
    C = [this.sp2.P[0]+.1,this.sp2.P[1],0];
    Q = [this.sp.P[0]-.1,this.sp.P[1]+.1,0];
    
    this.moveTentacle(2,C,R,Q,".51+(j+1)*.1");
    
    //animate mouth tentacles
    
    for(var i = 0 ; i < this.pTents.length ; i++){
        
        if(i<6){
            var string1 = "var fun2 = {y:.071,rz:-1*(.1+(this.mnoise(("+i+"+time*.6)+(-i*.051),0,.1)*(i/6))),rx:.1+(this.mnoise(("+i+"+1+time*.6)+(-i*.051),-.05,.05)*(i/6)),sx:.9,sy:.9,sz:.9}";
            var string2 = "var fun = {x:"+o.ball212.matrix[12]+"+.1,y:"+o.ball212.matrix[13]+",z:"+o.ball212.matrix[14]+",rz:3,ry:"+i/2+"-.06}";
            this.movePtentacle(this.pTents[i],this.spTents[i],string1,string2)
        }
        else{
            var string1 = "var fun2 = {y:.071,rz:"+msin(7*this.rotSpeed+(i*.5),-.02,.02)+"-.05,sx:.9,sy:.9,sz:.9,rot:true}";
            var string2 = "var fun = {x:"+o.ball25.matrix[12]+"+.1,y:"+o.ball25.matrix[13]+",z:"+o.ball25.matrix[14]+",rz:1.8,ry:"+i+"-3.14}";
            this.movePtentacle(this.pTents[i],this.spTents[i],string1,string2)
        }
    }
}

canvas1.movePtentacle=function(tent,sptent,func,func2){
    
    eval(func2);
    tent[0].recurse(makeMatrix(fun));
    
    var tempRZ = -1;
    
    for(var i = 0 ; i < tent.length ; i++){
        eval(func);
        var tempMat = makeMatrix(fun2);
        
        for (var j = i ; j > 0 ; j--){
            tempRZ-=(msin(7*this.rotSpeed+(i*.5),-.02,.02)-.05)/tent.length;
        }
        
        if(i>0){
            tent[i].recurse(tempMat);
            tent[i].func = fun2;
        }
        
        sptent[i].spring.F = [tent[i].matrix[12],tent[i].matrix[13],tent[i].matrix[14]];
        sptent[i].spring.update(map_range(i,0,tent.length,.1,.071));
        
        sptent[i].matrix = makeMatrix({
                                      rx:0,
                                      rz:tempRZ,
                                      x:sptent[i].spring.P[0],
                                      y:sptent[i].spring.P[1],
                                      z:sptent[i].spring.P[2]});
        
        // if(fun2.rot){
        //   sptent[i].matrix[0] = -tent[i].matrix[0];
        //   sptent[i].matrix[1] = -tent[i].matrix[1];
        //   sptent[i].matrix[2] = -tent[i].matrix[2];
        //   sptent[i].matrix[4] = -tent[i].matrix[4];
        //   sptent[i].matrix[5] = -tent[i].matrix[5];
        //   sptent[i].matrix[6] = -tent[i].matrix[6];
        //   sptent[i].matrix[8] = -tent[i].matrix[8];
        //   sptent[i].matrix[9] = -tent[i].matrix[9];
        //   sptent[i].matrix[10] = -tent[i].matrix[10];
        // }
    }
}

canvas1.makeUni=function(){
    
    this.unicycle = [];
    this.wheel = [];
    this.spokes = []
    this.pedals = [];
    
    //wheel parent
    
    this.addObject([0], 'fs_phong');
    var obj = this.objects[this.objects.length-1];
    this.things["wheel"] = obj;
    
    this.addObject([0], 'fs_phong');
    var obj = this.objects[this.objects.length-1];
    this.things["wangleParent"] = obj;
    
    this.things.wheel.add(this.things.wangleParent);
    
    //create outer tire
    var parms = {r1:1,r2:.1};
    this.addObject(createParametric(1/48,1/12,tor,parms), 'bike');
    var obj = this.objects[this.objects.length-1];
    obj.setUniform('p', [.02,.01, .01,0,0,0, 0,0,0,1]);
    obj.drawType="TRIANGLE_STRIP";
    this.things["outerTire"] = obj;
    this.wheel.push(obj);
    
    var parms = {r1:.95,r2:.09};
    this.addObject(createParametric(1/48,1/12,tor,parms), 'bike');
    var obj = this.objects[this.objects.length-1];
    obj.name = "innerTire";
    obj.setUniform('p', [.06,.04, .04,0,0,0, 0,0,0,1]);
    obj.drawType="TRIANGLE_STRIP";
    this.things["innerTire"] = obj;
    
    this.things.outerTire.add(obj);
    this.things.outerTire.add(this.things.innerTire);
    
    this.unicycle.push(this.things.outerTire);
    this.unicycle.push(this.things.innerTire);
    
    //make spokes
    for(var i = 0 ; i < 13; i ++){
        var parms = {r1:.02,r2:1};
        this.addObject(createParametric(1/9,1/7,cyl,parms), 'fs_phong');
        var obj = this.objects[this.objects.length-1];
        obj.setUniform('p', [.02,.01, .01,0,0,0, 0,0,0,1]);
        obj.drawType="TRIANGLE_STRIP";
        var name = "spoke"+i;
        this.things[name] = name;
        this.things.outerTire.add(obj);
        this.spokes.push(obj);
        
        
        this.unicycle.push(obj);
        
        
        for(var v = 0 ; v < obj.vertices.length ; v+=8){
            
            var rotparms =[];
            rotparms.rz = (i*Math.PI*2)/13;
            
            var rotMat = makeMatrix(rotparms);
            
            var vp = [];
            vp.x = obj.vertices[v];
            vp.y = obj.vertices[v+1];
            vp.z = obj.vertices[v+2];
            
            var vMat = makeMatrix(vp);
            
            var newMat = jMatMult(vMat,rotMat);
            
            obj.vertices[v+0] = newMat[12];
            obj.vertices[v+1] = newMat[13];
            obj.vertices[v+2] = newMat[14];
            
        }
        
        obj.vertexBuffer = createVertexBuffer(this.canvas.gl, obj.vertices);
        
    }
    
    //make axle
    
    this.makeThing(
                   {rx:Math.PI/2,sx:.1,sy:.1,sz:10,z:.05},
                   "axle",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:1,r2:.1}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.axle);
    
    
    this.unicycle.push(this.things.axle);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:.12,sy:.12,sz:1,z:.62,x:.02},
                   "smallAxle1",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:1,r2:.2}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.smallAxle1);
    this.unicycle.push(this.things.smallAxle1);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:.12,sy:.12,sz:1,z:-.62,x:.02,ry:Math.PI},
                   "smallAxle2",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:1,r2:.2}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.smallAxle2);
    this.unicycle.push(this.things.smallAxle2);
    
    //pedals
    
    this.makeThing(
                   {rz:Math.PI/2,sx:.12,sy:.12,sz:.12,z:4.4,x:.02,ry:Math.PI},
                   "jangle1",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:.6,r2:5}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.jangle1);
    this.unicycle.push(this.things.jangle1);
    
    
    this.makeThing(
                   {rz:Math.PI/2,sx:.12,sy:.12,sz:.12,z:-4.4,x:.02,ry:Math.PI*2},
                   "jangle2",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:.6,r2:5}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.jangle2);
    this.unicycle.push(this.things.jangle2);
    
    //pedal pivots
    this.makeThing(
                   {rx:-Math.PI/2,sx:.12,sy:.12,sz:.12,z:4,x:-5.02},
                   "jangle1A",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:.6,r2:5}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.jangle1A);
    this.unicycle.push(this.things.jangle1A);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:.12,sy:.12,sz:.12,z:-4,x:5.02},
                   "jangle2A",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   createParametric(1/24,1/18,cyl,{r1:.6,r2:5}),
                   'bike'
                   );
    this.things.outerTire.add(this.things.jangle2A);
    this.unicycle.push(this.things.jangle2A);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:1,sy:1,sz:1,x:-.3,z:.32},//,z:-4,x:5.02},
                   "pedalJoint1",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   [0],
                   'bike'
                   );
    this.things.outerTire.add(this.things.pedalJoint1);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:1,sy:1,sz:1,x:.3,z:-.32},//,z:-4,x:5.02},
                   "pedalJoint2",
                   [.06,.04, .04,0,0,0, 0,0,0,1],
                   [0],
                   'bike'
                   );
    this.things.outerTire.add(this.things.pedalJoint2);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:2,sy:.5,sz:3},//,x:-.3,z:.32},//,z:-4,x:5.02},
                   "pedal1",
                   [.5,.45, .49,0,0,0, 0,0,0,1],
                   createCube(),
                   'bike'
                   );
    this.things.pedalJoint1.add(this.things.pedal1);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:2,sy:.5,sz:3},//,x:-.3,z:.32},//,z:-4,x:5.02},
                   "pedal2",
                   [.5,.45, .49,0,0,0, 0,0,0,1],
                   createCube(),
                   'bike'
                   );
    this.things.pedalJoint2.add(this.things.pedal2);
    
    this.makeThing(
                   {rx:Math.PI/2,sx:.1,sy:.1,sz:2,z:.05},        //initial translation
                   "smallAxle",                                  //name
                   [.06,.04, .04,0,0,0, 0,0,0,1],                //uniform
                   createParametric(1/24,1/12,cyl,{r1:2,r2:.1}), //object
                   'bike'                                    //material
                   );
    this.things.outerTire.add(this.things.smallAxle);
    this.unicycle.push(this.things.smallAxle);
    
    this.makeThing(
                   {ry:Math.PI/2,rz:Math.PI/2,y:.45,sx:2,sy:2,sz:2},
                   "wangle",
                   [.5,.02, .02,0,0,0, 0,0,0,1],
                   createParametric(1/32,1/24,tor,{r1:.16,r2:.04}),
                   'bike'
                   );
    this.things.wangleParent.add(this.things.wangle);
    this.unicycle.push(this.things.wangle);
    
    this.makeThing(
                   {z:1},
                   "wangleton",
                   [.01,.08, .06,0,0,0, 0.02,0,0.025,1],
                   createParametric(1/92,1/14,tor,{r1:1.5,r2:.7}),
                   'fs_phong'
                   );
    
    this.makeThing(
                   {ry:Math.PI/2,rz:Math.PI/2,y:1.4,sx:2,sy:2,sz:2},
                   "seatPoint",
                   [1,.02, .02,0,0,0, 0,0,0,1],
                   [0],
                   'bike'
                   );
    this.things.wangleParent.add(this.things.seatPoint);
    
    //wingle the wangle
    var o = this.things.wangle;
    
    for(var i = 0 ; i < o.vertices.length ; i+=8){
        if(o.vertices[i]>0){
            o.vertices[i]+=5*.45;
            o.vertices[i]*=.2;
        }
    }
    
    o.vertexBuffer = createVertexBuffer(this.canvas.gl, o.vertices);
    this.things.wheel.add(this.things.outerTire);
    
    this.makeThing(
                   {y:.6,sx:2,sy:2,sz:2},
                   "shaft",
                   [.5,.02, .02,0,0,0, 0,0,0,1],
                   createParametric(1/12,1/18,cyl,{r1:.04,r2:.5}),
                   'bike'
                   );
    this.things.wangleParent.add(this.things.shaft);
    this.unicycle.push(this.things.shaft);
    
    for(t in this.unicycle){
        this.unicycle[t].textures = [];
        this.unicycle[t].textures[0] = "images/outdoor.jpg";
        this.unicycle[t].setUniform('bumpAmount',0.);
    }
    this.things.outerTire.setUniform('bumpAmount',1.);
}

canvas1.onmousemove = function(event) { // Mouse moved
    moveMouse(this.handle, event);
}

function moveMouse(handle, event) {
    var x = event.clientX;
    var y = event.clientY;
    var rect = event.target.getBoundingClientRect();
    if ( rect.left <= x && x <= rect.right &&
        rect.top  <= y && y <= rect.bottom ) {
        handle.mouseX = x - rect.left;
        handle.mouseY = y - rect.top;
    }
}

canvas1.moveTentacle = function(tent,C,R,Q,scale){
    
    if(this.tents[tent]){
        for(var j = 0 ; j < this.tents[tent].length ; j++){
            
            var obj = this.tents[tent][j];
            
            this.points = [];
            this.points.push(new Point(R));
            this.points.push(new Point(Q));
            this.points.push(new Point(C));
            
            var sp = new Quad(this.points);
            var parms = [];
            
            if(obj.name == ("ball"+tent+j)){
                var u = this.tents[tent].length;
                parms.y = sp.getPoint(j/u).y;
                parms.x = sp.getPoint(j/u).x;
                parms.z = sp.getPoint(j/u).z;
            }
            
            var scaleMat = makeMatrix(parms);
            obj.pMatrix = scaleMat;
            obj.matrix = recurse(obj,obj.pMatrix);
            
            var params = [];
            params.sx = params.sy = params.sz = eval(scale);
            var sMat = makeMatrix(params);
            
            obj.matrix = jMatMult(sMat,obj.matrix);
        }
    }
} 

canvas1.particles = function(num){
    
    var vee = this.things.particles.vertices;
    // console.log(this.things.mountains.vertices[this.things.mountains.vertices.length-7]);
    
    for(var i = 0 ; i <5 ; i++){
        this.things.particles.vertices.push(this.things.mountains.vertices[100]);
        this.things.particles.vertices.push(0);
        this.things.particles.vertices.push(-45);
        this.things.particles.vertices.push(1);
        this.things.particles.vertices.push(1);
        this.things.particles.vertices.push(1);
        this.things.particles.vertices.push(1);
        this.things.particles.vertices.push(1);
        this.prti.px.push(0);
        this.prti.py.push(0);
        this.prti.pz.push(0);
        this.prti.vx.push(0);
        this.prti.vy.push(0);
        this.prti.vz.push(0);
    }
    var b=0;
    for (var i = 0 ; i < this.things.particles.vertices.length ; i+=8){
        this.prti.px[b]  = map_range(this.n3d.noise3d( this.things.particles.vertices[i]+time*.1, this.things.particles.vertices[i+1], this.things.particles.vertices[i+2]),-1,1,0,.1);
        this.prti.px[b+1]=.1+this.n3d.noise3d( this.things.particles.vertices[i], this.things.particles.vertices[i+1], this.things.particles.vertices[i+2]);
        this.prti.pz[b+2]=.1*this.n3d.noise3d( this.things.particles.vertices[i], this.things.particles.vertices[i+1], this.things.particles.vertices[i+2]);
        
        this.prti.px[b]  = this.n3d.noise3d(time+b+i,0,0);
        this.prti.px[b+1]= .2;
        this.prti.pz[b+2]= 0;
        
        
        
        
        this.things.particles.vertices[i]  +=.1*this.prti.px[b] + (Math.random()*.01);
        this.things.particles.vertices[i+1]+=.1*this.prti.px[b+1];
        this.things.particles.vertices[i+2]+=.1*this.prti.px[b+2];
        
        b+=3;
        
        this.things.particles.vertices[i+3]-=Math.random()*.008;
        this.things.particles.vertices[i+4]-=Math.random()*.003;
        this.things.particles.vertices[i+5]-=Math.random()*.009;
    }
    this.things.particles.vertexBuffer = createVertexBuffer(this.canvas.gl, this.things.particles.vertices);
    while (this.things.particles.vertices.length>num){
        this.things.particles.vertices.shift();
    }
    /*
     this.prti = {
     px:[],
     py:[],
     pz:[],
     vx:[],
     vy:[],
     vz:[]
     };*/
}

window.onkeydown = onKeyDown;

function onKeyDown(evt) {
    console.log(evt.keyCode);
    //add a ball 'a'
    
    if(evt.keyCode == 81){
        canvas1.dance = !canvas1.dance;
        if(canvas1.dance)
            setSpan("foo","Dancing Jester");
        else
            setSpan("foo","Floppy Harlequin");
    }
    if(evt.keyCode == 69){
        canvas1.addSpeed += 1;
    }
    if(evt.keyCode == 68){
        d.debug=!d.debug;
    }
    if(evt.keyCode == 87){
        canvas1.addSpeed -= 1;
    }
    if(evt.keyCode == 68){
        canvas1.addSpeed *= 2;
    }
    if(evt.keyCode == 68){
        d.debug=!d.debug;
    }
    if(evt.keyCode == 83){
        canvas1.addSpeed /=2;
    }
}







