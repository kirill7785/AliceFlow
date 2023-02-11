unit UnitFlyuidMatLib;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, CheckLst;

type
  TFormFluidLibMat = class(TForm)
    GBFluidMatLib: TGroupBox;
    LFluidMatLib: TLabel;
    BApply: TButton;
    ComboBoxFluidLibMaterial: TComboBox;
    Label1: TLabel;
    Buttonconductivity: TButton;
    Label2: TLabel;
    Label3: TLabel;
    ButtonHeatCapacity: TButton;
    Label4: TLabel;
    Label5: TLabel;
    LabelDensityValue: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    ButtonDynamicViscosity: TButton;
    Label8: TLabel;
    Label9: TLabel;
    Buttonbeta: TButton;
    Label10: TLabel;
    procedure BApplyClick(Sender: TObject);
    procedure ButtonconductivityClick(Sender: TObject);
    procedure ButtonDynamicViscosityClick(Sender: TObject);
    procedure ButtonbetaClick(Sender: TObject);
    procedure ButtonHeatCapacityClick(Sender: TObject);
    procedure ComboBoxFluidLibMaterialChange(Sender: TObject);
    procedure FormClose(Sender: TObject; var Action: TCloseAction);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormFluidLibMat: TFormFluidLibMat;

implementation
 uses
     VisualUnit, addBlockUnit, Unitrectangularplot;
{$R *.dfm}

// Нажатие на кнопку Apply-применить
// Выбор билиотечного материала.
procedure TFormFluidLibMat.BApplyClick(Sender: TObject);
begin
    // 26.07.2016.
   // Выделен только один билиотечный материал
   with Laplas.body[Laplas.itek] do
   begin
      //
      Laplas.workmat[imatid].blibmat:=1; // означает что это библиотечный материал
      Laplas.workmat[imatid].ilibident:=1+ComboBoxFluidLibMaterial.ItemIndex; // номер библиотечного материала в соответствии с заданным
      Laplas.workmat[imatid].bBoussinesq:=1; // Для библиотечного материала приближение Обербека-Буссинеска  используется обязательным образом.
   end;
   Close;
end;

procedure TFormFluidLibMat.ButtonbetaClick(Sender: TObject);
begin
    frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;
   // linear expansion coefficient
    case ComboBoxFluidLibMaterial.ItemIndex of
   0 : begin
          // Dry Air
          frmRectangularPlot.cht1.Title.Text[0]:='Dry Air linear expansion coefficient';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='linear expansion coefficient x0.001, 1/K';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
          frmRectangularPlot.cht1.Series[0].AddXY(200,5.0325438,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,3.342765185,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,2.503025,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,2.002415,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,1.670511,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,1.434299,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
     1 : begin
          // Water Liquid
          frmRectangularPlot.cht1.Title.Text[0]:='Water Liquid linear expansion coefficient';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='linear expansion coefficient x0.0001, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
          frmRectangularPlot.cht1.Series[0].AddXY(273,-0.6875,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(278,0.09535,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(283,0.80397,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(293,2.028,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(313,3.8864,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(333,5.2476743,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(353,6.39073,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,7.5346,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
    end;
end;

procedure TFormFluidLibMat.ButtonconductivityClick(Sender: TObject);
begin
   frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;
   //  thermal conductivity
    case ComboBoxFluidLibMaterial.ItemIndex of
   0 : begin
          // Dry Air
          frmRectangularPlot.cht1.Title.Text[0]:='Dry Air conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity x0.01, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
          frmRectangularPlot.cht1.Series[0].AddXY(200,2.18364,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,2.522848,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,2.795017,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,3.02619,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,3.22919,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,3.411422,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
     1 : begin
          // Water Liquid
          frmRectangularPlot.cht1.Title.Text[0]:='Water Liquid conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity x0.1, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
          //frmRectangularPlot.cht1.Series[0].AddXY(700,1.7143,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(273,5.53,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(293,5.8618,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(313,6.1936,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(333,6.5254,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(353,6.8572,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,7.189,'',clRed);


          frmRectangularPlot.ShowModal;
       end;
    end;
end;

procedure TFormFluidLibMat.ButtonDynamicViscosityClick(Sender: TObject);
begin
    frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;
     //  dynamic viscosity
    case ComboBoxFluidLibMaterial.ItemIndex of
   0 : begin
          // Dry Air
          frmRectangularPlot.cht1.Title.Text[0]:='Dry Air dynamic viscosity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='dynamic viscosity x0.00001, Pa*s';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
          frmRectangularPlot.cht1.Series[0].AddXY(200,1.56536,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,1.765423,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,1.92268,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,2.05425,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,2.1684,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,2.26986,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
     1 : begin
          // Water Liquid
          frmRectangularPlot.cht1.Title.Text[0]:='Water Liquid dynamic viscosity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='dynamic viscosity x0.001, Pa*s';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
           frmRectangularPlot.cht1.Series[0].AddXY(273,1.801,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(293,1.012,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(313,0.65386,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(333,0.4583,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(353,0.339133,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,0.2609,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
    end;
end;

procedure TFormFluidLibMat.ButtonHeatCapacityClick(Sender: TObject);
begin
    frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;
    // Heat Capacity
     case ComboBoxFluidLibMaterial.ItemIndex of
   0 : begin
          // Dry Air
          frmRectangularPlot.cht1.Title.Text[0]:='Dry Air heat capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='heat capacity, J/(kgxK)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
          frmRectangularPlot.cht1.Series[0].AddXY(200,991.79223,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,1003.69622,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,1015.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,1027.504,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,1039.408,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,1051.31,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
     1 : begin
          // Water Liquid
          frmRectangularPlot.cht1.Title.Text[0]:='Water Liquid heat capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='heat capacity, J/(kgxK)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // соответствует базе данных AliceFlow_v0_24.
           frmRectangularPlot.cht1.Series[0].AddXY(273,4194,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(293,4177,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(313,4172,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(333,4179,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(353,4198,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,4229,'',clRed);

          frmRectangularPlot.ShowModal;
       end;
    end;
end;

procedure TFormFluidLibMat.ComboBoxFluidLibMaterialChange(Sender: TObject);
begin
   case ComboBoxFluidLibMaterial.ItemIndex of
   0 : begin
          // Dry Air
          LabelDensityValue.Caption:='1.10174';
       end;
   1 : begin
          // Water Liquid
          LabelDensityValue.Caption:='874.525';
       end;
   end;
end;

procedure TFormFluidLibMat.FormClose(Sender: TObject; var Action: TCloseAction);
begin
   // 26.07.2016.
   // Выделен только один билиотечный материал
   with Laplas.body[Laplas.itek] do
   begin
      //
      Laplas.workmat[imatid].blibmat:=1; // означает что это библиотечный материал
      Laplas.workmat[imatid].ilibident:=1+ComboBoxFluidLibMaterial.ItemIndex; // номер библиотечного материала в соответствии с заданным
      Laplas.workmat[imatid].bBoussinesq:=1; // Для библиотечного материала приближение Обербека-Буссинеска  используется обязательным образом.
   end;
end;

end.
