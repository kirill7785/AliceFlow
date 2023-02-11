unit addBlockUnit;
// редактирует добавляемый блок



interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TAddBlockForm = class(TForm)
    Panelglobalconteiner: TPanel;
    RadioGroupglobalconteiner: TRadioGroup;
    Panelinfo: TPanel;
    GBname: TGroupBox;
    Ename: TEdit;
    GroupBoxPriority: TGroupBox;
    EditPriority: TEdit;
    Bapply: TButton;
    btncolor: TButton;
    dlgColorcube: TColorDialog;
    lblset_trans: TLabel;
    lbltransparent: TLabel;
    scrlbrtrans1: TScrollBar;
    lblOpaque: TLabel;
    lblTransparentlab: TLabel;
    lbltransparencyvalue: TLabel;
    PanelGeometry: TPanel;
    GBsizeBlock: TGroupBox;
    LxS: TLabel;
    LyS: TLabel;
    LzS: TLabel;
    LxE: TLabel;
    LyE: TLabel;
    LzE: TLabel;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    ExS: TEdit;
    EyS: TEdit;
    EzS: TEdit;
    ExE: TEdit;
    EyE: TEdit;
    EzE: TEdit;
    PanelProperties: TPanel;
    RadioGroupType: TRadioGroup;
    Panel1: TPanel;
    Label7: TLabel;
    ButtonRadiation: TButton;
    GroupBoxPropBl: TGroupBox;
    LMN: TLabel;
    lblmateril: TLabel;
    GBmaterial: TGroupBox;
    RGSelect: TRadioGroup;
    BEditApply: TButton;
    CBselectAction: TComboBox;
    GBPower: TGroupBox;
    ButtonTransient: TButton;
    LabelGeometryType: TLabel;
    ComboBoxgeometrytype: TComboBox;
    LabelPlane: TLabel;
    ComboBoxPlane: TComboBox;
    Labelpowerinfo: TLabel;
    CheckBoxVisible: TCheckBox;
    GBPolygonGeom: TGroupBox;
    ListBoxvert: TListBox;
    ButtonAdd: TButton;
    ButtonRemove: TButton;
    EditHeight: TEdit;
    Editx: TEdit;
    Edity: TEdit;
    Editz: TEdit;
    Labeldim1: TLabel;
    Labeldimx: TLabel;
    Labeldimy: TLabel;
    Labeldimz: TLabel;
    LabelnameH: TLabel;
    Labelx: TLabel;
    Labely: TLabel;
    Labelz: TLabel;
    CheckBoxCylinder2Prism: TCheckBox;
    CheckBoxFixedCylinder: TCheckBox;
    LabelLineWidth: TLabel;
    ComboBoxLineWidth: TComboBox;
    GroupBoxActivityCondition: TGroupBox;
    ComboBoxZnak: TComboBox;
    LabelOperation: TLabel;
    LabelLeftExpression: TLabel;
    ComboBoxLeftExpression: TComboBox;
    EditLeftExpression: TEdit;
    EditRightExpression: TEdit;
    LabelRightExpression: TLabel;
    PanelNetwork: TPanel;
    GroupBoxinx: TGroupBox;
    ComboBoxinx: TComboBox;
    GroupBoxiny: TGroupBox;
    ComboBoxiny: TComboBox;
    GroupBoxinz: TGroupBox;
    ComboBoxinz: TComboBox;
    procedure BapplyClick(Sender: TObject);
    procedure BEditApplyClick(Sender: TObject);
    procedure RadioGroupTypeClick(Sender: TObject);
    procedure RGSelectClick(Sender: TObject);
    procedure btncolorClick(Sender: TObject);
    procedure scrlbrtrans1Change(Sender: TObject);
    procedure ButtonTransientClick(Sender: TObject);
    procedure ButtonRadiationClick(Sender: TObject);
    procedure RadioGroupglobalconteinerClick(Sender: TObject);
    procedure ComboBoxgeometrytypeChange(Sender: TObject);
    procedure ButtonRemoveClick(Sender: TObject);
    procedure ButtonAddClick(Sender: TObject);
    procedure ListBoxvertClick(Sender: TObject);
    procedure FormClose(Sender: TObject);
    procedure ComboBoxLeftExpressionChange(Sender: TObject);



    
  private
    { Private declarations }
     procedure setEdit(Sender: TObject);

  public
    { Public declarations }

  end;

var
  AddBlockForm: TAddBlockForm;

implementation
uses
     VisualUnit, UnitSolidMatLib, UnitUserDefinedSolidMaterial,
  UnitFlyuidMatLib, UnitUserDefinedFluidMaterial, UnitSelProjMat,
  UnitVariables, UnitTimedependpowerLaw, UnitBlockViewFactors;
{$R *.dfm}


// для полигона устанавливает значения значения окон редактирования.
procedure TAddBlockForm.setEdit(Sender: TObject);
begin
  if (ListBoxvert.ItemIndex<Laplas.body[Laplas.itek].nsizei) then
  begin
     EditHeight.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex]);
     Editx.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex]);
     Edity.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex]);
     Editz.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex]);
  end;
end;

// редактирование свойств добавляемого блока
procedure TAddBlockForm.BapplyClick(Sender: TObject);
var
   k,k1,i_1 : Integer; // текущий кубик
   // вспомогательные переменные для обработки исключительной ситуации
   bOk : Boolean;
   s1, s2, s3, s4, s5, s6,  sprior : String;
   r1, r2, r3, r4, r5, r6 : Real;
   priority_is_correct : Integer;

begin
    Laplas.bpolypoint:=false;

   // инициализация :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;
   r4:=0.0;
   r5:=0.0;
   r6:=0.0;

   // ввод данных
   k:=Laplas.itek;
   with Laplas.body[k] do
   begin
       BodyLineWidth:=1+ComboBoxLineWidth.ItemIndex; // толщина линии 1-6.


       bOk:=true; // признак правильности ввода

      // координаты блока

      if ((ComboBoxgeometrytype.ItemIndex=0)or(ComboBoxgeometrytype.ItemIndex=1)) then
      begin
         // Prism
      s1:=ExS.Text;  // параметризованные
      s2:=EyS.Text;  // геометрические
      s3:=ExE.Text;  // размеры
      s4:=EyE.Text;  // заданные
      s5:=EzS.Text;  // пользователем
      s6:=EzE.Text;

      if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
      if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
      if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учётом подстановки
      if bOk then r4:=FormVariables.my_real_convert(s4,bOk);  // значений переменных.
      if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
      if bOk then r6:=FormVariables.my_real_convert(s6,bOk);
      end;


      name:=Ename.Text; // имя блока


      // корректировка имени объекта чтобы избежать совпадающих имён.
      Laplas.correctobjname('b',name,k);
      Ename.Text:=name;
      AddBlockForm.Caption:='Blocks [ '+Trim(name)+' ]';

      sprior:= EditPriority.Text;
      for k1 := 1 to length(sprior) do
      begin
          if ((sprior[k1]<'0')or(sprior[k1]>'9')) then
          begin
            bOk:=false;
            EditPriority.Text:=IntToStr(priority);
            ShowMessage('Error! Priority is not positive integer value.');
          end;
      end;
      if (bOK) then
      begin
         priority_is_correct:=StrToInt(EditPriority.Text);
         if (bOk) then
         begin
            if (priority_is_correct>0) then
            begin
               // Приоритет может быть только положительным.
               // Жёсткое требование : приоритеты никаких двух блоков не могут совпадать.
               for k1 := 1 to Laplas.lb-1 do
               begin
                  if (k1<>k) then
                  begin
                      if (Laplas.body[k1].priority=priority_is_correct) then
                      begin
                         // Ошибка приоритеты двух существующих блоков совпали и не понятно кто
                         // кого перезаписывает при перекрытии.
                         bOk:=false;
                         EditPriority.Text:=IntToStr(priority);
                         ShowMessage('Error! compare priority in blocks:'+Laplas.body[k1].name+' and '+name);
                      end;
                  end;
              end;
            end
             else
            begin
               // Приоритет не может быть отрицательным.
               // Присваиваем полю ввода реальный текущий приоритет блока.
               EditPriority.Text:=IntToStr(priority);
               bOk:=false;
               ShowMessage('Error! Priority must be positive value.');
            end;
          end;
      end;




      if (bOk) then
      begin



         priority:=priority_is_correct;
         if (priority>Laplas.priority_id) then
         begin
            // Если введённый пользователем приоритет
            // выше любого существующего в программе на данный момент
            // то мы корректируем наибольшее значение приоритта которое будет присвоено новому блоку.
            Laplas.priority_id:=priority+5;
         end;

         if (ComboBoxgeometrytype.ItemIndex=0) then
         begin
            // Prism
            CylinderFixed:=false;

            igeometry_type:=0;
            iPlane:=1; // XY

            sxS:=s1;  // параметризованные
            syS:=s2;  // геометрические
            sxE:=s3;  // размеры
            syE:=s4;  // заданные
            szS:=s5;  // пользователем
            szE:=s6;


            xS:=r1;  // числовые размеры
            yS:=r2;  // заданные пользователем
            xE:=r3;  // с учётом подстановки
            yE:=r4;  // значений переменных.
            zS:=r5;
            zE:=r6;

            xC:=0.5*(abs(xE+xS));
            yC:=0.5*(abs(yE+yS));
            zC:=zS;
            Hcyl:=abs(zE-zS);
            R_out_cyl:=0.5*(0.5*(abs(xE-xS))+0.5*(abs(yE-yS)));
            R_in_cyl:=0.0;


            sxC:=FormatFloat('0.000',xC);
            syC:=FormatFloat('0.000',yC);
            szC:=FormatFloat('0.000',zC);
            sHcyl:=FormatFloat('0.000',Hcyl);
            sR_out_cyl:=FormatFloat('0.000',R_out_cyl);
            sR_in_cyl:=FormatFloat('0.000',R_in_cyl);

            // Нужно ли преобразовывать цилиндр в призму.
            bCylinder2Prism:=CheckBoxCylinder2Prism.Checked;

         end;

         if (ComboBoxgeometrytype.ItemIndex=2) then
         begin
            // Polygon
            CylinderFixed:=false;

            igeometry_type:=2;
            iPlane_obj2:=ComboBoxPlane.ItemIndex+1;
            s1:=EditHeight.Text;  // параметризованные
            s2:=Editx.Text;  // геометрические
            s3:=Edity.Text;  // размеры
            s4:=Editz.Text;

            bOk:=true;
            if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
            if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
            if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учётом подстановки
            if bOk then r4:=FormVariables.my_real_convert(s4,bOk);  // значений переменных.

            if (bOk) then
            begin
              xi[ListBoxvert.ItemIndex]:=r2;
              yi[ListBoxvert.ItemIndex]:=r3;
              zi[ListBoxvert.ItemIndex]:=r4;
              hi[ListBoxvert.ItemIndex]:=r1;
            end;
            for i_1 := 0 to ListBoxvert.Count-1 do
            begin
               hi[i_1]:=r1;
            end;

            case iPlane_obj2 of
              1 : begin
                     // XY
                     for i_1 := 0 to ListBoxvert.Count-1 do
                     begin
                        zi[i_1]:=r4;
                     end;
                end;
                 2 : begin
                     // XZ
                     for i_1 := 0 to ListBoxvert.Count-1 do
                     begin
                        yi[i_1]:=r3;
                     end;
                end;
                3 : begin
                     // YZ
                     for i_1 := 0 to ListBoxvert.Count-1 do
                     begin
                        xi[i_1]:=r2;
                     end;
                end;
            end;

            // Нужно ли преобразовывать цилиндр в призму.
            bCylinder2Prism:=CheckBoxCylinder2Prism.Checked;

         end;

         // Cylindr
         if (ComboBoxgeometrytype.ItemIndex=1) then
         begin

            if (CheckBoxFixedCylinder.Checked) then
            begin
               CylinderFixed:=true;
            end
            else
            begin
               CylinderFixed:=false;
            end;

            // Cylinder
            igeometry_type:=1;
            iPlane:=ComboBoxPlane.ItemIndex+1; // плоскость основания цилиндра выбранная пользователем.

            sxC:=s1;  // параметризованные
            syC:=s2;  // геометрические
            szC:=s5;  // размеры
            sHcyl:=s3;  // заданные
            sR_out_cyl:=s4;  // пользователем
            sR_in_cyl:=s6;


            xC:=r1;  // числовые размеры
            yC:=r2;  // заданные пользователем
            zC:=r5;  // с учётом подстановки
            Hcyl:=r3;  // значений переменных.
            R_out_cyl:=r4;
            R_in_cyl:=r6;

            case iPlane of
               1 : begin
                      // XY
                      xS:=xC-R_out_cyl;
                      yS:=yC-R_out_cyl;
                      xE:=xC+R_out_cyl;
                      yE:=yC+R_out_cyl;
                      zS:=zC;
                      zE:=zC+Hcyl;
                   end;
               2 : begin
                      // XZ
                      xS:=xC-R_out_cyl;
                      zS:=zC-R_out_cyl;
                      xE:=xC+R_out_cyl;
                      zE:=zC+R_out_cyl;
                      yS:=yC;
                      yE:=yC+Hcyl;
                   end;
               3 : begin
                      // YZ
                      yS:=yC-R_out_cyl;
                      zS:=zC-R_out_cyl;
                      yE:=yC+R_out_cyl;
                      zE:=zC+R_out_cyl;
                      xS:=xC;
                      xE:=xC+Hcyl;
                   end;
            end;

            sxS:=FormatFloat('0.000',xS);
            syS:=FormatFloat('0.000',yS);
            szS:=FormatFloat('0.000',zS);
            sxE:=FormatFloat('0.000',xE);
            syE:=FormatFloat('0.000',yE);
            szE:=FormatFloat('0.000',zE);

            // Нужно ли преобразовывать цилиндр в призму.
            bCylinder2Prism:=CheckBoxCylinder2Prism.Checked;

         end;

         if (ComboBoxgeometrytype.ItemIndex=3) then
         begin
            // 18.09.2022
            // CAD
            CylinderFixed:=false;

            igeometry_type:=3;
         end;

          // Activity condition.
            // 0 ==
            // 1 !=
            // 2 <=
            // 3 <
            // 4 >=
            // 5 >
            iOperation:=ComboBoxZnak.ItemIndex;
            // 0 expression
            // 1 X Min
            // 2 X Max
            // 3 Y Min
            // 4 Y Max
            // 5 Z Min
            // 6 Z Max
            iLeftExpression:=ComboBoxLeftExpression.ItemIndex;
            sLeftOperation:=Trim(EditLeftExpression.Text);
            sRightOperation:=Trim(EditRightExpression.Text);

            if ((igeometry_type=0)or(igeometry_type=2)) then
            begin
                // PRISM
                // Или POLYGON
                inx_network_conteiner:=ComboBoxinx.ItemIndex+1;
                iny_network_conteiner:=ComboBoxiny.ItemIndex+1;
                inz_network_conteiner:=ComboBoxinz.ItemIndex+1;
            end
            else
            begin
               // Геометрия пользователя отличная от кубика.
               inx_network_conteiner:=1;
               iny_network_conteiner:=1;
               inz_network_conteiner:=1;
            end;


      end;

      transparency:=1.0-1.0*(scrlbrtrans1.Position-scrlbrtrans1.Min)/(scrlbrtrans1.Max-scrlbrtrans1.Min);
      bvisible:=CheckBoxVisible.Checked;
   end;

   if (bOk) then
   begin
      Laplas.body[Laplas.itek].bactivity:=Laplas.activity_body(Laplas.body[Laplas.itek]);
   end;

   with Laplas do
   begin
      ReadyPaint;
   end;
end;

// Редактирует или назначает свойства материала.
procedure TAddBlockForm.BEditApplyClick(Sender: TObject);
var
    i : Integer;
    iold_matid : Integer;
begin
   if ((CBselectAction.ItemIndex=0) or (CBselectAction.ItemIndex=1)) then
   begin
      // 0 - Edit Current Material
      if (CBselectAction.ItemIndex=1) then
      begin
         FormUserDefinedSolidMat.ButtonCancel.Visible:=true;
         FormUserDefinedFluidMaterial.ButtonCancel.Visible:=true;

         // 1 - Create New Material
         if ((RadioGroupType.ItemIndex=0) or (RadioGroupType.ItemIndex=2)) then
         begin
            iold_matid:=Laplas.body[Laplas.itek].imatid;
            Laplas.body[Laplas.itek].imatid:=Laplas.lmatmax;
            inc(Laplas.lmatmax);
            SetLength(Laplas.workmat,Laplas.lmatmax);
            with Laplas.workmat[Laplas.lmatmax-1] do
            begin
               if (RadioGroupType.ItemIndex=0) then
               begin


                  //SOLID
                  rho:=2800; // плотность дюр-аллюминия
                  //cp:=921; // удельная теплоёмкость дюр-аллюминия
                  n_cp:=1;
                  SetLength(temp_cp, n_cp);
                  SetLength(arr_cp, n_cp);
                  temp_cp[0]:=20.0;
                  arr_cp[0]:= 921; // теплоёмкость дюр-аллюминия
                  //lambda:=164; // теплопроводность дюр-аллюминия
                  n_lam:=1;
                  SetLength(temp_lam, n_lam);
                  SetLength(arr_lam, n_lam);
                  temp_lam[0]:=20.0;
                  arr_lam[0]:= 164; // теплопроводность дюр-аллюминия

                  // Ортотропность теплопроводности.
                  mult_lam_x:=1.0;
                  mult_lam_y:=1.0;
                  mult_lam_z:=1.0;

                  // Stress
                  n_Poisson_ratio:=1;
                  SetLength(temp_Poisson_ratio,n_Poisson_ratio);
                  SetLength(arr_Poisson_ratio,n_Poisson_ratio);
                  temp_Poisson_ratio[0]:=20.0;
                  arr_Poisson_ratio[0]:=0.334;
                  n_Young_Module:=1;
                  SetLength(temp_Young_Module,n_Young_Module);
                  SetLength(arr_Young_Module,n_Young_Module);
                  temp_Young_Module[0]:=20.0;
                  arr_Young_Module[0]:=69.0;  // GPa aluminium
                  n_Linear_expansion_coefficient:=1;
                  SetLength(temp_Linear_expansion_coefficient, n_Linear_expansion_coefficient);
                  SetLength(arr_Linear_expansion_coefficient, n_Linear_expansion_coefficient);
                  temp_Linear_expansion_coefficient[0]:=20.0;
                  arr_Linear_expansion_coefficient[0]:=23.0;

                  mult_Linear_expansion_coefficient_x:=1.0;
                  mult_Linear_expansion_coefficient_y:=1.0;
                  mult_Linear_expansion_coefficient_z:=1.0;
                  mult_Young_Module_x:=1.0;
                  mult_Young_Module_y:=1.0;
                  mult_Young_Module_z:=1.0;
                  mult_Poisson_ratio_xy:=1.0;
                  mult_Poisson_ratio_xz:=1.0;
                  mult_Poisson_ratio_yz:=1.0;
                  mult_Poisson_ratio_yx:=1.0;
                  mult_Poisson_ratio_zx:=1.0;
                  mult_Poisson_ratio_zy:=1.0;
                  bShearModuleActive:=false;
                  FormUserDefinedSolidMat.CheckBoxShearModulus.Checked:=false;
                  FormUserDefinedSolidMat.CheckBoxShearModulusClick(Sender);
                  ShearModuleGxy:=arr_Young_Module[0]/(2.0*(1.0+arr_Poisson_ratio[0]));
                  ShearModuleGyz:=arr_Young_Module[0]/(2.0*(1.0+arr_Poisson_ratio[0]));
                  ShearModuleGxz:=arr_Young_Module[0]/(2.0*(1.0+arr_Poisson_ratio[0]));


                  // следующие два значения не используются
                  mu:=1.7894e-5; // воздух
                  beta_t:=0.003331; // воздух
                  bBoussinesq:=0; // приближение Обербека-Буссинеска выключено.
                  namemat:='Al-Duralumin';
                  blibmat:=0; // материал определяемый пользователем.
                  ilibident:=100; // в библиотеке материалов не значится.
                  ilawmu:=0; // const
                  mumin:=mu; mumax:=mu; // ограничители вязкости
                  Amu:=1.0; Bmu:=1.0; Cmu:=1.0; // параметры модели для вязкости
                  degreennmu:=1.0; // показатель степени
               end;
               if (RadioGroupType.ItemIndex=2) then
               begin
                  //FLUID
                  rho:=1.1614; // плотность воздуха
                  //cp:=1005; // удельная теплоёмкость воздуха
                  n_cp:=1;
                  SetLength(temp_cp, n_cp);
                  SetLength(arr_cp, n_cp);
                  temp_cp[0]:=20.0;
                  arr_cp[0]:= 1005; // теплоёмкость воздуха
                  // lambda:=0.0261; // теплопроводность воздуха
                  n_lam:=1;
                  SetLength(temp_lam, n_lam);
                  SetLength(arr_lam, n_lam);
                  temp_lam[0]:=20.0;
                  arr_lam[0]:= 0.0261; // теплопроводность воздуха
                  // Ортотропность теплопроводности.
                  mult_lam_x:=1.0;
                  mult_lam_y:=1.0;
                  mult_lam_z:=1.0;

                  // Stress
                  n_Poisson_ratio:=1;
                  SetLength(temp_Poisson_ratio,n_Poisson_ratio);
                  SetLength(arr_Poisson_ratio,n_Poisson_ratio);
                  temp_Poisson_ratio[0]:=20.0;
                  arr_Poisson_ratio[0]:=0.49;  // air
                  n_Young_Module:=1;
                  SetLength(temp_Young_Module,n_Young_Module);
                  SetLength(arr_Young_Module,n_Young_Module);
                  temp_Young_Module[0]:=20.0;
                  arr_Young_Module[0]:=1.42e-4; // GPa air
                  n_Linear_expansion_coefficient:=1;
                  SetLength(temp_Linear_expansion_coefficient, n_Linear_expansion_coefficient);
                  SetLength(arr_Linear_expansion_coefficient, n_Linear_expansion_coefficient);
                  temp_Linear_expansion_coefficient[0]:=20.0;
                  arr_Linear_expansion_coefficient[0]:=3.331e+3; // air

                  mult_Linear_expansion_coefficient_x:=1.0;
                  mult_Linear_expansion_coefficient_y:=1.0;
                  mult_Linear_expansion_coefficient_z:=1.0;
                  mult_Young_Module_x:=1.0;
                  mult_Young_Module_y:=1.0;
                  mult_Young_Module_z:=1.0;
                  mult_Poisson_ratio_xy:=1.0;
                  mult_Poisson_ratio_xz:=1.0;
                  mult_Poisson_ratio_yz:=1.0;
                  mult_Poisson_ratio_yx:=1.0;
                  mult_Poisson_ratio_zx:=1.0;
                  mult_Poisson_ratio_zy:=1.0;
                  bShearModuleActive:=false;
                  FormUserDefinedSolidMat.CheckBoxShearModulus.Checked:=false;
                  FormUserDefinedSolidMat.CheckBoxShearModulusClick(Sender);
                  ShearModuleGxy:=arr_Young_Module[0]/(2.0*(1.0+arr_Poisson_ratio[0]));
                  ShearModuleGyz:=arr_Young_Module[0]/(2.0*(1.0+arr_Poisson_ratio[0]));
                  ShearModuleGxz:=arr_Young_Module[0]/(2.0*(1.0+arr_Poisson_ratio[0]));


                  mu:=1.7894e-5; // коэффициент динамической вязкости
                  beta_t:=0.003331; // коэффициент линейного температурного расширения
                  namemat:='air';
                  bBoussinesq:=0; // приближение Обербека-Буссинеска выключено
                  blibmat:=0; // это материал определённый пользователем
                  ilibident:=0; // это материал определённый пользователем
                  ilawmu:=0; // const
                  mumin:=mu; mumax:=mu; // ограничители вязкости
                  Amu:=1.0; Bmu:=1.0; Cmu:=1.0;  // параметры модели для вязкости
                  degreennmu:=1.0; // показатель степени
               end;
            end;
         end;
      end
      else
      begin
        // Отмена может быть только при создании нового материала.
        FormUserDefinedSolidMat.ButtonCancel.Visible:=false;
        FormUserDefinedFluidMaterial.ButtonCancel.Visible:=false;
      end;


   if (RadioGroupType.ItemIndex=0) then
   begin
      //SOLID
      if (RGSelect.ItemIndex=0) then
      begin
         if (CBselectAction.ItemIndex=1) then
         begin
            // Create new matherial
             FormSolidLibMat.cbbSolidMatLib.ItemIndex:=0;
             FormSolidLibMat.cbbSolidMatLib.OnChange(Sender);
             // Мы создаём новый библиотечный твёрдотельный материал которого раньше не было.
             // Новые настройки всегда изотропны.
             FormSolidLibMat.Editx.Text:='1';
             FormSolidLibMat.Edity.Text:='1';
             FormSolidLibMat.Editz.Text:='1';
         end;

         // SOLID Program Library
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=1) then
         begin
            // Это вызывается при редактировании существующего библиотечного матеиала.
            // инициализация
            FormSolidLibMat.cbbSolidMatLib.ItemIndex:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilibident-101;
            FormSolidLibMat.cbbSolidMatLib.OnChange(Sender);
            // Ортотропность теплопроводности.
            FormSolidLibMat.Editx.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_lam_x);
            FormSolidLibMat.Edity.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_lam_y);
            FormSolidLibMat.Editz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_lam_z);
         end;
         if ((CBselectAction.ItemIndex=0)and((Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=0))) then
         begin
           // Текущий не библиотечный материал нельзя редактировать как библиотечный.
           RGSelect.ItemIndex:=1; // переключаем с библиотечного на пользовательский
           // материал.
         end
         else
         begin
            // Истинно библиотечный материал.
         FormSolidLibMat.ShowModal; // библиотечные материалы
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=1) then
         begin
            case Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilibident of
              101 : begin
                    LMN.Caption:='Alumina';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Alumina';
                    end;
              102 : begin
                    LMN.Caption:='Si';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Si';
                    end;
              103 : begin
                    LMN.Caption:='GaAs';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='GaAs';
                    end;
              104 : begin
                    LMN.Caption:='GaN';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='GaN';
                    end;
              105 : begin
                    LMN.Caption:='SiC4H';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='SiC4H';
                    end;
              106 : begin
                    LMN.Caption:='Sapphire';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Sapphire';
                    end;
              107 : begin
                    LMN.Caption:='Diamond';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Diamond';
                    end;
              108 : begin
                    LMN.Caption:='MD40';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='MD40';
                    end;
              109 : begin
                    LMN.Caption:='Au';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Au';
                    end;
              110 : begin
                    LMN.Caption:='SiO2';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='SiO2';
                    end;
              111 : begin
                    LMN.Caption:='Cu';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Cu';
                    end;
              112 : begin
                    LMN.Caption:='Kovar';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Kovar';
                    end;
              113 : begin
                    LMN.Caption:='Brass LS 59-1-L';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Brass_LS_59_1_L';
                    end;
              114 : begin
                    LMN.Caption:='Al-Duralumin';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Al_Duralumin';
                    end;
              115 : begin
                    LMN.Caption:='AlN';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='AlN';
                    end;
            end;
         end;
         end;
      end;
      // User-Defined material
      if (RGSelect.ItemIndex=1) then
      begin
         FormUserDefinedSolidMat.CBsolidmat.ItemIndex:=0; // no pattern
         // SOLID User - Defined
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=0) then
         begin
            // инициализация:
            FormUserDefinedSolidMat.EMatName.Text:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
            //FormUserDefinedSolidMat.EditPoissonRatio.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].Poisson_ratio);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_Poisson_ratio=1) then
            begin
               // Constant properties.
               FormUserDefinedSolidMat.ComboBoxPoissonratio.ItemIndex:=0;
               FormUserDefinedSolidMat.ComboBoxPoissonratioClick(Sender);
               FormUserDefinedSolidMat.EditPoissonRatio.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_Poisson_ratio[0]);
            end
             else
            begin
               // Piecewise Poissonratio properties.
               FormUserDefinedSolidMat.ComboBoxPoissonratio.ItemIndex:=1;
               FormUserDefinedSolidMat.ComboBoxPoissonratioClick(Sender);
            end;
            //FormUserDefinedSolidMat.EditYoungModule.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].Young_Module);
             if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_Young_Module=1) then
            begin
               // Constant properties.
               FormUserDefinedSolidMat.ComboBoxYoungModule.ItemIndex:=0;
               FormUserDefinedSolidMat.ComboBoxYoungModuleClick(Sender);
               FormUserDefinedSolidMat.EditYoungModule.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_Young_Module[0]);
            end
             else
            begin
               // Piecewise EditYoungModule properties.
               FormUserDefinedSolidMat.ComboBoxYoungModule.ItemIndex:=1;
               FormUserDefinedSolidMat.ComboBoxYoungModuleClick(Sender);
            end;
            //FormUserDefinedSolidMat.EditLinearExpansionKoefficient.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].Linear_expansion_coefficient);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_Linear_expansion_coefficient=1) then
            begin
               // Constant properties.
               FormUserDefinedSolidMat.ComboBoxlinearExpansion.ItemIndex:=0;
               FormUserDefinedSolidMat.ComboBoxlinearExpansionChange(Sender);
               FormUserDefinedSolidMat.EditLinearExpansionKoefficient.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_Linear_expansion_coefficient[0]);
            end
             else
            begin
               // Piecewise linear properties.
               FormUserDefinedSolidMat.ComboBoxlinearExpansion.ItemIndex:=1;
               FormUserDefinedSolidMat.ComboBoxlinearExpansionChange(Sender);
            end;
            FormUserDefinedSolidMat.Erho.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].rho);
            //FormUserDefinedSolidMat.ECp.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].cp);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp=1) then
            begin
               // Constant properties.
               FormUserDefinedSolidMat.ComboBoxheatcapacitytype.ItemIndex:=0;
               FormUserDefinedSolidMat.ComboBoxheatcapacitytypeChange(Sender);
               FormUserDefinedSolidMat.ECp.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[0]);
            end
             else
            begin
               // Piecewise linear properties.
               FormUserDefinedSolidMat.ComboBoxheatcapacitytype.ItemIndex:=1;
               FormUserDefinedSolidMat.ComboBoxheatcapacitytypeChange(Sender);
            end;

            //FormUserDefinedSolidMat.ELam.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].lambda);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam=1) then
            begin
               // Constant properties.
               FormUserDefinedSolidMat.ComboBoxconductivitytype.ItemIndex:=0;
               FormUserDefinedSolidMat.ComboBoxconductivitytypeChange(Sender);
               FormUserDefinedSolidMat.ELam.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[0]);
            end
             else
            begin
               // Piecewise linear properties.
               FormUserDefinedSolidMat.ComboBoxconductivitytype.ItemIndex:=1;
               FormUserDefinedSolidMat.ComboBoxconductivitytypeChange(Sender);
            end;

            FormUserDefinedSolidMat.Editmultx.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_lam_x);
            FormUserDefinedSolidMat.Editmulty.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_lam_y);
            FormUserDefinedSolidMat.Editmultz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_lam_z);

            FormUserDefinedSolidMat.EditbetaX.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Linear_expansion_coefficient_x);
            FormUserDefinedSolidMat.EditbetaY.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Linear_expansion_coefficient_y);
            FormUserDefinedSolidMat.EditbetaZ.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Linear_expansion_coefficient_z);

            FormUserDefinedSolidMat.EditEx.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Young_Module_x);
            FormUserDefinedSolidMat.EditEy.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Young_Module_y);
            FormUserDefinedSolidMat.EditEz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Young_Module_z);

            FormUserDefinedSolidMat.Editnuxy.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Poisson_ratio_xy);
            FormUserDefinedSolidMat.Editnuxz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Poisson_ratio_xz);
            FormUserDefinedSolidMat.Editnuyz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mult_Poisson_ratio_yz);

            
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].bShearModuleActive) then
            begin
               FormUserDefinedSolidMat.CheckBoxShearModulus.Checked:=true;
            end
             else
            begin
               FormUserDefinedSolidMat.CheckBoxShearModulus.Checked:=false;
            end;

            FormUserDefinedSolidMat.CheckBoxShearModulusClick(Sender);

            FormUserDefinedSolidMat.EditGxy.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].ShearModuleGxy);
            FormUserDefinedSolidMat.EditGyz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].ShearModuleGyz);
            FormUserDefinedSolidMat.EditGxz.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].ShearModuleGxz);

         end;
         if (Laplas.bREALESEversion) then
         begin
            // МАНИФЕСТ СТАБИЛЬНОСТИ
            // ДЛЯ ПОЛЬЗОВАТЕЛЯ 04.08.2019

            // Мы полностью отключаем все
            // недоработанные функциональные
            // возможности из интерфейса программы.


            // Не задаём свойства механических материалов.
            //FormUserDefinedSolidMat.GroupBoxThermalStress.Visible:=false;
            //FormUserDefinedSolidMat.PanelSOLID.Width:=287;
            //FormUserDefinedSolidMat.Width:=290;

         end;

         FormUserDefinedSolidMat.bCanselSelect:=false;
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].lambda_is_file=1) then
         begin
            FormUserDefinedSolidMat.CheckBox_lambda_is_file.Checked:=true;
            FormUserDefinedSolidMat.GBuserproperties.Visible:=false;
         end
          else
         begin
            FormUserDefinedSolidMat.CheckBox_lambda_is_file.Checked:=false;
            FormUserDefinedSolidMat.GBuserproperties.Visible:=true;
         end;

          if (Laplas.egddata.iStaticStructural=1) then
          begin
             // Механика активна.
             FormUserDefinedSolidMat.ClientWidth:=614;
             FormUserDefinedSolidMat.GroupBoxThermalStress.Visible:=true;
          end
           else
          begin
             // Механика неактивна.
             FormUserDefinedSolidMat.ClientWidth:=287;
             FormUserDefinedSolidMat.GroupBoxThermalStress.Visible:=false;
          end;

         FormUserDefinedSolidMat.ShowModal;
         if ((CBselectAction.ItemIndex=1) and (FormUserDefinedSolidMat.bCanselSelect)) then
         begin
            // была отмена  создания нового материала.
            Laplas.body[Laplas.itek].imatid:=iold_matid;
            LMN.Caption:=Laplas.workmat[iold_matid].namemat;
            // Новый материал не будет создан.
            dec(Laplas.lmatmax);
            SetLength(Laplas.workmat,Laplas.lmatmax);
         end
         else
         begin
            LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
         end;
      end;
   end;
   if (RadioGroupType.ItemIndex=2) then
   begin
      //FLUID
      if (RGSelect.ItemIndex=0) then
      begin
         // FLUID Program Library
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=1) then
         begin
            // инициализация
            FormFluidLibMat.ComboBoxFluidLibMaterial.ItemIndex:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilibident-1;
            case FormFluidLibMat.ComboBoxFluidLibMaterial.ItemIndex of
              0 : begin
                     // Dry Air
                     FormFluidLibMat.LabelDensityValue.Caption:='1.10174';
                  end;
              1 : begin
                     // Water Liquid
                     FormFluidLibMat.LabelDensityValue.Caption:='874.525';
                  end;
            end;
         end;
         FormFluidLibMat.ShowModal; // библиотечные материалы
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=1) then
         begin
            case Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilibident of
              1 : begin
                    LMN.Caption:='Dry_Air';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Dry_Air';
                    end;
              2 : begin
                    LMN.Caption:='Water_Liquid';
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Water_Liquid';
                    end;
              end;
         end;
      end;
      // User-Defined Fluid material
      if (RGSelect.ItemIndex=1) then
      begin
         FormUserDefinedFluidMaterial.CBTipPattern.ItemIndex:=0; // no pattern
          // FLUID User - Defined
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=0) then
         begin
            // инициализация:
            FormUserDefinedFluidMaterial.CBRho.ItemIndex:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].bBoussinesq;
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].bBoussinesq=1) then
            begin
               // делаем видимым задание линейного температурного коэффициента расширения.
               FormUserDefinedFluidMaterial.GBVolexpans.Visible:=true;
            end
            else
            begin
               // делаем невидимым задание линейного температурного коэффициента расширения.
               FormUserDefinedFluidMaterial.GBVolexpans.Visible:=false;
            end;

            FormUserDefinedFluidMaterial.EMatName.Text:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
            FormUserDefinedFluidMaterial.ERho.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].rho);
            //FormUserDefinedFluidMaterial.ECp.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].cp);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp=1) then
            begin
               // Constant properties
               FormUserDefinedFluidMaterial.CBCp.ItemIndex:=0;
               FormUserDefinedFluidMaterial.CBCpChange(Sender);
               FormUserDefinedFluidMaterial.ECp.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[0]);
            end
             else
            begin
               // Piecewise linear properties.
               FormUserDefinedFluidMaterial.CBCp.ItemIndex:=1;
               FormUserDefinedFluidMaterial.CBCpChange(Sender);
            end;
            //FormUserDefinedFluidMaterial.ELam.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].lambda);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam=1) then
            begin
               // Constant properties.
               FormUserDefinedFluidMaterial.CBLam.ItemIndex:=0;
               FormUserDefinedFluidMaterial.CBLamChange(Sender);
               FormUserDefinedFluidMaterial.ELam.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[0]);
            end
             else
            begin
              // piecewise properties
              FormUserDefinedFluidMaterial.CBLam.ItemIndex:=1;
              FormUserDefinedFluidMaterial.CBLamChange(Sender);
            end;
            // динамическая вязкость
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilawmu=0) then
            begin
               FormUserDefinedFluidMaterial.CBMu.ItemIndex:=0;
               FormUserDefinedFluidMaterial.LSIMu.Visible:=true;
               FormUserDefinedFluidMaterial.EMu.Visible:=true;
               FormUserDefinedFluidMaterial.EMu.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mu);
            end
            else
            begin
               FormUserDefinedFluidMaterial.CBMu.ItemIndex:=1;
               FormUserDefinedFluidMaterial.EMu.Visible:=false;
               FormUserDefinedFluidMaterial.LSIMu.Visible:=false;
               FormUserDefinedFluidMaterial.EMu.Text:=' ';
            end;

            FormUserDefinedFluidMaterial.EBeta_T.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].beta_t);
         end;
         FormUserDefinedFluidMaterial.bCanselSelect:=false;
         FormUserDefinedFluidMaterial.ShowModal;
         if ((CBselectAction.ItemIndex=1) and (FormUserDefinedFluidMaterial.bCanselSelect)) then
         begin
            // была отмена  создания нового материала.
            Laplas.body[Laplas.itek].imatid:=iold_matid;
            LMN.Caption:=Laplas.workmat[iold_matid].namemat;
            // Новый материал не будет создан.
            dec(Laplas.lmatmax);
            SetLength(Laplas.workmat,Laplas.lmatmax);
         end
         else
         begin
            LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
         end;
      end;
   end;
   end;

   if (CBselectAction.ItemIndex=2) then
   begin
      // Select Project Material
      FormSelProjMat.CBlistprojmat.Clear; // очистка списка
      for i:=0 to Laplas.lmatmax-1 do
      begin
         // заполнение списка
         FormSelProjMat.CBlistprojmat.Items.Append(Laplas.workmat[i].namemat);
      end;
      FormSelProjMat.CBlistprojmat.ItemIndex:=0;
      FormSelProjMat.ShowModal; // Вызов формы
      case Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat of
      0 : begin
             // материал заданный пользователем
             RGSelect.ItemIndex:=1;
          end;
      1 : begin
             // библиотечный материал
             RGSelect.ItemIndex:=0;
          end;
      end; // case
      LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat; // печатаем имя выбранного материала
   end;
end;


procedure TAddBlockForm.RadioGroupglobalconteinerClick(Sender: TObject);
begin
   case RadioGroupglobalconteiner.ItemIndex of
      0 : begin
             // Info
             PanelGeometry.Visible:=false;
             Panelinfo.Visible:=true;
             PanelProperties.Visible:=false;
             PanelNetwork.Visible:=false;
      end;
      1 : begin
             // Geometry
             PanelGeometry.Visible:=true;
             Panelinfo.Visible:=false;
             PanelProperties.Visible:=false;
             PanelNetwork.Visible:=false;
      end;
      2 : begin
             // Properties
             PanelGeometry.Visible:=false;
             Panelinfo.Visible:=false;
             PanelProperties.Visible:=true;
             PanelNetwork.Visible:=false;
             if (ComboBoxgeometrytype.ItemIndex=2) then
             begin
                // Polygon
                Laplas.body[Laplas.itek].n_power:=1;
                SetLength(Laplas.body[Laplas.itek].arr_s_power, Laplas.body[Laplas.itek].n_power);
                SetLength(Laplas.body[Laplas.itek].arr_power, Laplas.body[Laplas.itek].n_power);
                SetLength(Laplas.body[Laplas.itek].temp_power, Laplas.body[Laplas.itek].n_power);
                if (FormatSettings.DecimalSeparator=',') then
                begin
                   Laplas.body[Laplas.itek].arr_s_power[0]:='0,0';
                end;
                if (FormatSettings.DecimalSeparator='.') then
                begin
                   Laplas.body[Laplas.itek].arr_s_power[0]:='0.0';
                end;
                Laplas.body[Laplas.itek].arr_power[0]:=0.0;
                Laplas.body[Laplas.itek].temp_power[0]:=20.0;
                if (FormatSettings.DecimalSeparator=',') then
                begin
                   LabelPowerInfo.Caption:='0,0';
                end;
                if (FormatSettings.DecimalSeparator='.') then
                begin
                   LabelPowerInfo.Caption:='0.0';
                end;
             end;
      end;
      3 : begin
             // Network
             PanelGeometry.Visible:=false;
             Panelinfo.Visible:=false;
             PanelProperties.Visible:=false;
             if (RadioGroupType.ItemIndex=0) then
             begin
                PanelNetwork.Visible:=true;
             end
             else
             begin
                PanelNetwork.Visible:=false;
             end;
        end;
   end;
end;

// Смена типа материала : SOLID, FLUID or Hollow.
procedure TAddBlockForm.RadioGroupTypeClick(Sender: TObject);
begin
   case RadioGroupType.ItemIndex of
    0 : begin
           // SOLID
            ButtonRadiation.Visible:=true;
           if (Laplas.body[Laplas.itek].itype=2) then
           begin
              // бывший Hollow Block (материала может не быть.
              Laplas.MainMemo.Lines.Add('материал '+IntToStr(Laplas.body[Laplas.itek].imatid));

           end;
           Laplas.body[Laplas.itek].itype:=1;

           if (Laplas.body[Laplas.itek].imatid<Length(Laplas.workmat)) then
           begin
              LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
           end
           else
           begin
              Laplas.body[Laplas.itek].imatid:=0;
              LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
           end;
           GroupBoxPropBl.Visible:=true;
        end;
    1 : begin
           // HOLLOW
           ButtonRadiation.Visible:=false;
           // нет внутренней surface - 2- surface  модели.
           Laplas.body[Laplas.itek].binternalRadiation:=0;
           Laplas.body[Laplas.itek].itype:=2;
           GroupBoxPropBl.Visible:=false;
        end;
    2 : begin
           // FLUID
           ButtonRadiation.Visible:=true;
           if (Laplas.body[Laplas.itek].itype=2) then
           begin
              // бывший Hollow Block (материала может не быть.
              Laplas.MainMemo.Lines.Add('материал '+IntToStr(Laplas.body[Laplas.itek].imatid));
           end;
           Laplas.body[Laplas.itek].itype:=3;
           if (Laplas.body[Laplas.itek].imatid<Length(Laplas.workmat)) then
           begin
              LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
           end
           else
           begin
              Laplas.body[Laplas.itek].imatid:=0;
              LMN.Caption:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
           end;
           GroupBoxPropBl.Visible:=true;
        end;
    end;
end;

// Реакция на смену схемы определения свойств материала :
// библиотека или определение пользователем.
procedure TAddBlockForm.RGSelectClick(Sender: TObject);
begin
    // Лучше не печатать это сообщение т.к. оно надоедает
   //Application.MessageBox('Change happens only when you click Edit','Attantion!')
end;


// Задаёт цвет кубу.
procedure TAddBlockForm.btncolorClick(Sender: TObject);
begin
   if dlgColorcube.Execute then
   begin
      Laplas.body[Laplas.itek].dcol:=dlgColorcube.Color;
      Laplas.ColorToGL(dlgColorcube.Color,Laplas.body[Laplas.itek].dcol, Laplas.body[Laplas.itek].redcolor,Laplas.body[Laplas.itek].greencolor,Laplas.body[Laplas.itek].bluecolor);
   end;
   Laplas.body[Laplas.itek].transparency:=1.0-1.0*(scrlbrtrans1.Position-scrlbrtrans1.Min)/(scrlbrtrans1.Max-scrlbrtrans1.Min);
end;

// Добавление точки в полигон.
procedure TAddBlockForm.ButtonAddClick(Sender: TObject);
var
   i : Integer;
begin
     SetLength(Laplas.body[Laplas.itek].xi,ListBoxvert.Count+1);
      SetLength(Laplas.body[Laplas.itek].yi,ListBoxvert.Count+1);
      SetLength(Laplas.body[Laplas.itek].zi,ListBoxvert.Count+1);
      SetLength(Laplas.body[Laplas.itek].hi,ListBoxvert.Count+1);
      for i := ListBoxvert.Count  to ListBoxvert.ItemIndex+1 do
      begin
         Laplas.body[Laplas.itek].xi[i]:=Laplas.body[Laplas.itek].xi[i-1];
         Laplas.body[Laplas.itek].yi[i]:=Laplas.body[Laplas.itek].yi[i-1];
         Laplas.body[Laplas.itek].zi[i]:=Laplas.body[Laplas.itek].zi[i-1];
         Laplas.body[Laplas.itek].hi[i]:=Laplas.body[Laplas.itek].hi[i-1];
      end;
      if (ListBoxvert.ItemIndex>0) then
      begin
         if (ListBoxvert.ItemIndex<ListBoxvert.Count-1) then
         begin
            Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex-1]+Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex+1]);
            Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex-1]+Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex+1]);
            Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex-1]+Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex+1]);
            Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex-1]+Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex+1]);
         end
         else
         begin
            Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex-1] + Laplas.body[Laplas.itek].xi[0]);
            Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex-1] + Laplas.body[Laplas.itek].yi[0]);
            Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex-1] + Laplas.body[Laplas.itek].zi[0]);
            Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex-1] + Laplas.body[Laplas.itek].hi[0]);
         end;
      end
      else
      begin
         Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex+1] + Laplas.body[Laplas.itek].xi[ListBoxvert.Count-1]);
         Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex+1] + Laplas.body[Laplas.itek].yi[ListBoxvert.Count-1]);
         Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex+1] + Laplas.body[Laplas.itek].zi[ListBoxvert.Count-1]);
         Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex]:=0.5*(Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex+1] + Laplas.body[Laplas.itek].hi[ListBoxvert.Count-1]);
      end;
      ListBoxvert.Items.Add('vert '+IntToStr(ListBoxvert.Count+1));
      inc(Laplas.body[Laplas.itek].nsizei);
      setEdit(Sender);
      // Загрузить новые координаты объекта.
end;

// Задаём зависимость мощности тепловыделения от времени.
procedure TAddBlockForm.ButtonRadiationClick(Sender: TObject);
begin
   // Показывает View Factors.
   // инициализация.
   //FormRadiation.Edit1.Text:=FloatToStr(Laplas.body[Laplas.itek].emissW);
   //FormRadiation.Edit2.Text:=FloatToStr(Laplas.body[Laplas.itek].emissE);
   //FormRadiation.Edit3.Text:=FloatToStr(Laplas.body[Laplas.itek].emissS);
   //FormRadiation.Edit4.Text:=FloatToStr(Laplas.body[Laplas.itek].emissN);
   //FormRadiation.Edit5.Text:=FloatToStr(Laplas.body[Laplas.itek].emissB);
   //FormRadiation.Edit6.Text:=FloatToStr(Laplas.body[Laplas.itek].emissT);

   FormRadiation.Edit1.Text:=Laplas.body[Laplas.itek].semissW;
   FormRadiation.Edit2.Text:=Laplas.body[Laplas.itek].semissE;
   FormRadiation.Edit3.Text:=Laplas.body[Laplas.itek].semissS;
   FormRadiation.Edit4.Text:=Laplas.body[Laplas.itek].semissN;
   FormRadiation.Edit5.Text:=Laplas.body[Laplas.itek].semissB;
   FormRadiation.Edit6.Text:=Laplas.body[Laplas.itek].semissT;

   if (Laplas.body[Laplas.itek].binternalRadiation=0) then
   begin
      FormRadiation.CheckBoxinternalRadiation.Checked:=false;
   end
    else
   begin
      FormRadiation.CheckBoxinternalRadiation.Checked:=true;
   end;
   
   FormRadiation.ShowModal;
end;

// Удаление опорной точки из полигона.
procedure TAddBlockForm.ButtonRemoveClick(Sender: TObject);
var
   i, id : Integer;
begin
   // Только если количество опорных точек больше трёх.
   if (Laplas.body[Laplas.itek].nsizei>3) then
   begin
      if (ListBoxvert.ItemIndex=ListBoxvert.Count-1) then
      begin
         id:=ListBoxvert.Count-2;
      end
      else
      begin
         id:=ListBoxvert.ItemIndex;
      end;
      for i := ListBoxvert.ItemIndex to ListBoxvert.Count-2 do
      begin
         Laplas.body[Laplas.itek].xi[i]:=Laplas.body[Laplas.itek].xi[i+1];
         Laplas.body[Laplas.itek].yi[i]:=Laplas.body[Laplas.itek].yi[i+1];
         Laplas.body[Laplas.itek].zi[i]:=Laplas.body[Laplas.itek].zi[i+1];
         Laplas.body[Laplas.itek].hi[i]:=Laplas.body[Laplas.itek].hi[i+1];
      end;
      SetLength(Laplas.body[Laplas.itek].xi,ListBoxvert.Count-1);
      SetLength(Laplas.body[Laplas.itek].yi,ListBoxvert.Count-1);
      SetLength(Laplas.body[Laplas.itek].zi,ListBoxvert.Count-1);
      SetLength(Laplas.body[Laplas.itek].hi,ListBoxvert.Count-1);
      ListBoxvert.Items.Delete(ListBoxvert.ItemIndex);
      dec(Laplas.body[Laplas.itek].nsizei);
      // Загрузить новые координаты объекта.
      ListBoxvert.ItemIndex:=id;
      setEdit(Sender);
   end;
end;

procedure TAddBlockForm.ButtonTransientClick(Sender: TObject);
begin
    // 2 - это полигон. Мы не можем задаваать мощность тепловыделения
    // внутри полигона.
    // 24.10.2020 Мощность внутри полигона можно задавать.
    //if (Laplas.body[Laplas.itek].igeometry_type<>2) then
    begin
       if (Laplas.body[Laplas.itek].n_power=1) then
       begin
          // Constant (not temperature depend)
          FormTransientPowerSetting.EditPower.Visible:=true;
          FormTransientPowerSetting.Label1.Visible:=true;
          FormTransientPowerSetting.Buttonpiecewisepower.Visible:=false;
          //FormTransientPowerSetting.EditPower.Text:=FloatToStr(Laplas.body[Laplas.itek].arr_power[0]);
          FormTransientPowerSetting.EditPower.Text:=Trim(Laplas.body[Laplas.itek].arr_s_power[0]);
          FormTransientPowerSetting.ComboBoxTemperaturedependpower.ItemIndex:=0;
          FormTransientPowerSetting.ComboBoxTemperaturedependpowerChange(Sender);
       end
        else
       begin
          // Temperature depend.
          FormTransientPowerSetting.EditPower.Visible:=false;
          FormTransientPowerSetting.Label1.Visible:=false;
          FormTransientPowerSetting.Buttonpiecewisepower.Visible:=true;
          FormTransientPowerSetting.ComboBoxTemperaturedependpower.ItemIndex:=1;
          FormTransientPowerSetting.ComboBoxTemperaturedependpowerChange(Sender);
       end;
       // Зависимость мощности тепловыделения от времени.
       // 0 - не зависит от времени и выделяется постоянно,
       // 1 - squre wave зависимость от времени,
       // 2 - square wave 2 зависимость от времени,
       // 3 - hot cold режим для Евдокимовой Н.Л.
       FormTransientPowerSetting.RadioGroupTimeDependPowerLow.ItemIndex:=Laplas.body[Laplas.itek].ipower_time_depend;
       FormTransientPowerSetting.ShowModal;
       //labelpowerinfo.Caption:=FloatToStr(Laplas.body[Laplas.itek].arr_power[0])+' W';
    end
    //else
    //begin
      // ShowMessage('Can not power define in polygon. 26.01.2018');
      // Laplas.MainMemo.Lines.Add('Can not power define in polygon. 26.01.2018');
    //end;
end;

// Смена типа геометрии.
procedure TAddBlockForm.ComboBoxgeometrytypeChange(Sender: TObject);
var
    i : Integer;
begin
   // Смена типа геометрии.
   case ComboBoxgeometrytype.ItemIndex of
      0 : begin
             // Prism
             GBPolygonGeom.Visible:=false;
             GBsizeBlock.Visible:=true;
             CheckBoxFixedCylinder.Visible:=false;
             CheckBoxFixedCylinder.Checked:=false;
             Laplas.body[Laplas.itek].CylinderFixed:=false;

             LxS.Caption:='xS';
             LyS.Caption:='yS';
             LzS.Caption:='zS';
             LxE.Caption:='xE';
             LyE.Caption:='yE';
             LzE.Caption:='zE';
             LabelPlane.Visible:=false;
             ComboBoxPlane.Visible:=false;

             // Может быть здесь более уместна параметризованная информация в виде стрингов на основе переменных.
             ExS.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].xS);  // параметризованные
             EyS.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].yS);  // геометрические
             ExE.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].xE);  // размеры
             EyE.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].yE);  // заданные
             EzS.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].zS);  // пользователем
             EzE.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].zE);

          end;
      1 : begin
             // Cylinder
             GBPolygonGeom.Visible:=false;
             GBsizeBlock.Visible:=true;

             if (Laplas.bREALESEversion) then
             begin
                // МАНИФЕСТ СТАБИЛЬНОСТИ
                // ДЛЯ ПОЛЬЗОВАТЕЛЯ 04.08.2019

                // Мы полностью отключаем все
                // недоработанные функциональные
                // возможности из интерфейса программы.
                // Убираем управление фиксацией боковых стенок цилиндра.
                CheckBoxFixedCylinder.Visible:=false;
             end
              else
             begin
                CheckBoxFixedCylinder.Visible:=true;
             end;
             Laplas.body[Laplas.itek].CylinderFixed:=false;
             CheckBoxFixedCylinder.Checked:=false;

             LxS.Caption:='xC';
             LyS.Caption:='yC';
             LzS.Caption:='zC';
             LxE.Caption:='Height';
             LyE.Caption:='Radius';
             LzE.Caption:='Int radius';
             LabelPlane.Visible:=true;
             ComboBoxPlane.Visible:=true;
             ComboBoxPlane.ItemIndex:=Laplas.body[Laplas.itek].iPlane-1; // 1 - XY

             // Может быть здесь более уместна параметризованная информация в виде стрингов на основе переменных.
              ExS.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].xC);  // параметризованные
              EyS.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].yC);  // геометрические
              ExE.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].Hcyl);  // размеры
              EyE.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].R_out_cyl);  // заданные
              EzS.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].zC);  // пользователем
              EzE.Text:=FormatFloat('0.000',Laplas.body[Laplas.itek].R_in_cyl);
          end;
      2 : begin
             // Polygon
             GBPolygonGeom.Visible:=true;
             GBsizeBlock.Visible:=false;
             CheckBoxFixedCylinder.Visible:=false;
             Laplas.body[Laplas.itek].CylinderFixed:=false;
             CheckBoxFixedCylinder.Checked:=false;

             LabelPlane.Visible:=true;
             ComboBoxPlane.Visible:=true;
             ComboBoxPlane.ItemIndex:=Laplas.body[Laplas.itek].iPlane_obj2-1; // 1 - XY

             ListBoxvert.Items.Clear;
             for i := 1 to Laplas.body[Laplas.itek].nsizei do
             begin
                ListBoxvert.Items.Add('vert '+IntToStr(i));
             end;
             ListBoxvert.ItemIndex:=0;
             setEdit(Sender);
             Labelx.Caption:='x1';
             Labely.Caption:='y1';
             Labelz.Caption:='z1';

             if (Laplas.ComboBoxlength.ItemIndex=1) then
             begin
                Labeldim1.Caption:='mm';
                Labeldimx.Caption:='mm';
                Labeldimy.Caption:='mm';
                Labeldimz.Caption:='mm';
             end;
             if (Laplas.ComboBoxlength.ItemIndex=0) then
             begin
                Labeldim1.Caption:='m';
                Labeldimx.Caption:='m';
                Labeldimy.Caption:='m';
                Labeldimz.Caption:='m';
             end;
             if (Laplas.ComboBoxlength.ItemIndex=2) then
             begin
                Labeldim1.Caption:='um';
                Labeldimx.Caption:='um';
                Labeldimy.Caption:='um';
                Labeldimz.Caption:='um';
             end;
          end;
          3 :
          begin
             // CAD
             LabelPlane.Visible:=false;
             ComboBoxPlane.Visible:=false;
             GBPolygonGeom.Visible:=false;
             GBsizeBlock.Visible:=false;
             CheckBoxFixedCylinder.Visible:=false;
             CheckBoxFixedCylinder.Checked:=false;
             Laplas.body[Laplas.itek].CylinderFixed:=false;
          end;
   end;
end;





procedure TAddBlockForm.ComboBoxLeftExpressionChange(Sender: TObject);
begin
   if (ComboBoxLeftExpression.ItemIndex=0) then
   begin
      EditLeftExpression.Visible:=true;
   end
   else
   begin
      EditLeftExpression.Visible:=false;
   end;
end;

procedure TAddBlockForm.FormClose(Sender: TObject);
begin
    Laplas.bpolypoint:=false;
    BapplyClick(Sender);
end;

// Загружает данные о координатах в поля для редактирования.
procedure TAddBlockForm.ListBoxvertClick(Sender: TObject);
begin
      EditHeight.Text:=FloatToStr(Laplas.body[Laplas.itek].hi[ListBoxvert.ItemIndex]);  // параметризованные
      Editx.Text:=FloatToStr(Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex]);  // геометрические
      Edity.Text:=FloatToStr(Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex]);  // размеры
      Editz.Text:=FloatToStr(Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex]);
      Laplas.xpolypoint:=Laplas.body[Laplas.itek].xi[ListBoxvert.ItemIndex];
      Laplas.ypolypoint:=Laplas.body[Laplas.itek].yi[ListBoxvert.ItemIndex];
      Laplas.zpolypoint:=Laplas.body[Laplas.itek].zi[ListBoxvert.ItemIndex];
      Laplas.bpolypoint:=true;
end;

procedure TAddBlockForm.scrlbrtrans1Change(Sender: TObject);
begin
   lbltransparencyvalue.Caption:=FloatToStr(1.0*(scrlbrtrans1.Position-scrlbrtrans1.Min)/(scrlbrtrans1.Max-scrlbrtrans1.Min));
   Laplas.body[Laplas.itek].transparency:=1.0-1.0*(scrlbrtrans1.Position-scrlbrtrans1.Min)/(scrlbrtrans1.Max-scrlbrtrans1.Min);
   with Laplas do
   begin
      ReadyPaint;
   end;
end;

end.
